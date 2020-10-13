#!/usr/bin/env julia --project
#=
# This experiment file establishes the initial conditions, boundary conditions,
# source terms and simulation parameters (domain size + resolution) for the
# SOARES LES case. The set of parameters presented in the `master` branch copy
# include those that have passed offline tests at the full simulation time of
# 6 hours. Suggested offline tests included plotting horizontal-domain averages
# of key properties (see AtmosDiagnostics). The timestepper configuration is in
# `src/Driver/solver_configs.jl` while the `AtmosModel` defaults can be found in
# `src/Atmos/Model/AtmosModel.jl` and `src/Driver/driver_configs.jl`
#
# This setup works in both Float32 and Float64 precision. `FT`
#
# To simulate the full 6 hour experiment, change `timeend` to (3600*6) and type in
#
# julia --project experiments/AtmosLES/soares.jl
#
# See `src/Driver/driver_configs.jl` for additional flags (e.g. VTK, diagnostics,
# update-interval, output directory settings)
#
# Upcoming changes:
# 1) Atomic sources
# 2) Improved boundary conditions
# 3) Collapsed experiment design
# 4) Updates to generally keep this in sync with master

@article{doi:10.1175/1520-0469(2003)60<1201:ALESIS>2.0.CO;2,
author = {Siebesma, A. Pier and Bretherton,
          Christopher S. and Brown,
          Andrew and Chlond,
          Andreas and Cuxart,
          Joan and Duynkerke,
          Peter G. and Jiang,
          Hongli and Khairoutdinov,
          Marat and Lewellen,
          David and Moeng,
          Chin-Hoh and Sanchez,
          Enrique and Stevens,
          Bjorn and Stevens,
          David E.},
title = {A Large Eddy Simulation Intercomparison Study of Shallow Cumulus Convection},
journal = {Journal of the Atmospheric Sciences},
volume = {60},
number = {10},
pages = {1201-1219},
year = {2003},
doi = {10.1175/1520-0469(2003)60<1201:ALESIS>2.0.CO;2},
URL = {https://journals.ametsoc.org/doi/abs/10.1175/1520-0469%282003%2960%3C1201%3AALESIS%3E2.0.CO%3B2},
eprint = {https://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%282003%2960%3C1201%3AALESIS%3E2.0.CO%3B2}
=#

using ArgParse
using Distributions
using DocStringExtensions
using LinearAlgebra
using Printf
using Random
using StaticArrays
using Test

using ClimateMachine
using ClimateMachine.Atmos
using ClimateMachine.Orientations
using ClimateMachine.ConfigTypes
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.ODESolvers
using ClimateMachine.Thermodynamics
using ClimateMachine.TurbulenceClosures
using ClimateMachine.TurbulenceConvection
using ClimateMachine.VariableTemplates
using ClimateMachine.BalanceLaws:
    BalanceLaw, Auxiliary, Gradient, GradientFlux, Prognostic

using CLIMAParameters
using CLIMAParameters.Planet: e_int_v0, grav, day
using CLIMAParameters.Atmos.Microphysics

struct LiquidParameterSet <: AbstractLiquidParameterSet end
struct IceParameterSet <: AbstractIceParameterSet end

struct MicropysicsParameterSet{L, I} <: AbstractMicrophysicsParameterSet
    liq::L
    ice::I
end

struct EarthParameterSet{M} <: AbstractEarthParameterSet
    microphys::M
end

microphys = MicropysicsParameterSet(LiquidParameterSet(), IceParameterSet())

const param_set = EarthParameterSet(microphys)

import ClimateMachine.Atmos: atmos_source!
using ClimateMachine.Atmos: altitude, recover_thermo_state

"""
  Soares Geostrophic Forcing (Source)
"""
struct SoaresGeostrophic{FT} <: Source
    "Coriolis parameter [s⁻¹]"
    f_coriolis::FT
    "Eastward geostrophic velocity `[m/s]` (Base)"
    u_geostrophic::FT
    "Eastward geostrophic velocity `[m/s]` (Slope)"
    u_slope::FT
    "Northward geostrophic velocity `[m/s]`"
    v_geostrophic::FT
end
function atmos_source!(
    s::SoaresGeostrophic,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)

    f_coriolis = s.f_coriolis
    u_geostrophic = s.u_geostrophic
    u_slope = s.u_slope
    v_geostrophic = s.v_geostrophic

    z = altitude(atmos, aux)
    # Note z dependence of eastward geostrophic velocity
    u_geo = SVector(u_geostrophic + u_slope * z, v_geostrophic, 0)
    ẑ = vertical_unit_vector(atmos, aux)
    fkvector = f_coriolis * ẑ
    # Accumulate sources
    source.ρu -= fkvector × (state.ρu .- state.ρ * u_geo)
    return nothing
end

"""
  Soares Sponge (Source)
"""
struct SoaresSponge{FT} <: Source
    "Maximum domain altitude (m)"
    z_max::FT
    "Altitude at with sponge starts (m)"
    z_sponge::FT
    "Sponge Strength 0 ⩽ α_max ⩽ 1"
    α_max::FT
    "Sponge exponent"
    γ::FT
    "Eastward geostrophic velocity `[m/s]` (Base)"
    u_geostrophic::FT
    "Eastward geostrophic velocity `[m/s]` (Slope)"
    u_slope::FT
    "Northward geostrophic velocity `[m/s]`"
    v_geostrophic::FT
end
function atmos_source!(
    s::SoaresSponge,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)

    z_max = s.z_max
    z_sponge = s.z_sponge
    α_max = s.α_max
    γ = s.γ
    u_geostrophic = s.u_geostrophic
    u_slope = s.u_slope
    v_geostrophic = s.v_geostrophic

    z = altitude(atmos, aux)
    u_geo = SVector(u_geostrophic + u_slope * z, v_geostrophic, 0)
    ẑ = vertical_unit_vector(atmos, aux)
    # Accumulate sources
    if z_sponge <= z
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α_max * sinpi(r / 2)^s.γ
        source.ρu -= β_sponge * (state.ρu .- state.ρ * u_geo)
    end
    return nothing
end

"""
  SoaresTendencies (Source)
No tendencies are applied in Soares
"""

"""
  Initial Condition for Soares LES
"""
function init_soares!(problem, bl, state, aux, (x, y, z), t)
    # This experiment runs the SOARES LES Configuration
    # (Shallow cumulus cloud regime)
    # x,y,z imply eastward, northward and altitude coordinates in `[m]`

    # Problem floating point precision
    FT = eltype(state)

    P_sfc::FT = 1.0e5 # Surface air pressure
    qg::FT = 5.0e-3 # Total moisture at surface
    q_pt_sfc = PhasePartition(qg) # Surface moisture partitioning
    Rm_sfc = gas_constant_air(bl.param_set, q_pt_sfc) # Moist gas constant
    θ_liq_sfc = FT(300) # Prescribed θ_liq at surface
    T_sfc = FT(300) # Surface temperature
    _grav = FT(grav(bl.param_set))

    # Initialise speeds [u = Eastward, v = Northward, w = Vertical]
    u::FT = 0.01
    v::FT = 0
    w::FT = 0

    # Assign piecewise quantities to θ_liq and q_tot
    θ_liq::FT = 0
    q_tot::FT = 0

    # Piecewise functions for potential temperature and total moisture
    if FT(0) <= z <= FT(1350)
        # Well mixed layer
        θ_liq = FT(300)
        q_tot = 5.0e-3 - 3.7e-4* (z / FT(1000))
    else
        # Conditionally unstable layer
        θ_liq = 300.0 + 2.0 * (z-1350.0)/1000.0
        q_tot = 5.0e-3 - 3.7e-4 * 1.35 - 9.4e-4 * (z-1350.0)/1000.0
    end
    # Convert total specific humidity to kg/kg
    q_tot /= 1000
    # Scale height based on surface parameters
    H = Rm_sfc * T_sfc / _grav
    # Pressure based on scale height
    P = P_sfc * exp(-z / H)

    # Establish thermodynamic state and moist phase partitioning
    TS = LiquidIcePotTempSHumEquil_given_pressure(bl.param_set, θ_liq, P, q_tot)
    T = air_temperature(TS)
    ρ = air_density(TS)
    q_pt = PhasePartition(TS)

    # Compute momentum contributions
    ρu = ρ * u
    ρv = ρ * v
    ρw = ρ * w

    # Compute energy contributions
    e_kin = FT(1 // 2) * (u^2 + v^2 + w^2)
    e_pot = _grav * z
    ρe_tot = ρ * total_energy(e_kin, e_pot, TS)

    # Assign initial conditions for prognostic state variables
    state.ρ = ρ
    state.ρu = SVector(ρu, ρv, ρw)
    state.ρe = ρe_tot
    state.moisture.ρq_tot = ρ * q_tot
    if bl.moisture isa NonEquilMoist
        state.moisture.ρq_liq = FT(0)
        state.moisture.ρq_ice = FT(0)
    end

    if z <= FT(400) # Add random perturbations to bottom 400m of model
        state.ρe += rand() * ρe_tot / 100
        state.moisture.ρq_tot += rand() * ρ * q_tot / 100
    end
    init_state_prognostic!(bl.turbconv, bl, state, aux, (x, y, z), t)
end

function soares_model(
    ::Type{FT},
    config_type,
    zmax,
    surface_flux;
    turbconv = NoTurbConv(),
    moisture_model = "equilibrium",
) where {FT}

    ics = init_soares!     # Initial conditions

    C_smag = FT(0.23)     # Smagorinsky coefficient

    u_star = FT(0.28)     # Friction velocity
    C_drag = FT(0.0011)   # Bulk transfer coefficient

    T_sfc = FT(300)     # Surface temperature `[K]`
    q_sfc = FT(5.0e-3)  # Surface specific humiity `[kg/kg]`
    LHF = FT(70)       # Latent heat flux `[W/m²]`
    SHF = FT(70)         # Sensible heat flux `[W/m²]`
    moisture_flux = LHF / latent_heat_vapor(param_set, T_sfc)

    z_sponge = FT(2400)     # Start of sponge layer
    α_max = FT(0.75)        # Strength of sponge layer (timescale)
    γ = 2              # Strength of sponge layer (exponent)

    u_geostrophic = FT(0.01)        # Eastward relaxation speed
    u_slope = FT(0)     # Slope of altitude-dependent relaxation speed
    v_geostrophic = FT(0)          # Northward relaxation speed
    f_coriolis = FT(0.376e-4) # Coriolis parameter

    # Assemble source components
    source_default = (
        Gravity(),
        SoaresSponge{FT}(
            zmax,
            z_sponge,
            α_max,
            γ,
            u_geostrophic,
            u_slope,
            v_geostrophic,
        ),
        SoaresGeostrophic{FT}(f_coriolis, u_geostrophic, u_slope, v_geostrophic),
        turbconv_sources(turbconv)...,
    )
    if moisture_model == "equilibrium"
        source = source_default
        moisture = EquilMoist{FT}(; maxiter = 5, tolerance = FT(0.1))
    elseif moisture_model == "nonequilibrium"
        source = (source_default..., CreateClouds())
        moisture = NonEquilMoist()
    else
        @warn @sprintf(
            """
%s: unrecognized moisture_model in source terms, using the defaults""",
            moisture_model,
        )
        source = source_default
    end

    # Set up problem initial and boundary conditions
    if surface_flux == "prescribed"
        energy_bc = PrescribedEnergyFlux((state, aux, t) -> LHF + SHF)
        moisture_bc = PrescribedMoistureFlux((state, aux, t) -> moisture_flux)
    elseif surface_flux == "bulk"
        energy_bc = BulkFormulaEnergy(
            (state, aux, t, normPu_int) -> C_drag,
            (state, aux, t) -> (T_sfc, q_sfc),
        )
        moisture_bc = BulkFormulaMoisture(
            (state, aux, t, normPu_int) -> C_drag,
            (state, aux, t) -> q_sfc,
        )
    else
        @warn @sprintf(
            """
%s: unrecognized surface flux; using 'prescribed'""",
            surface_flux,
        )
    end

    problem = AtmosProblem(
        boundarycondition = (
            AtmosBC(
                momentum = Impenetrable(DragLaw(
                    # normPu_int is the internal horizontal speed
                    # P represents the projection onto the horizontal
                    (state, aux, t, normPu_int) -> (u_star / normPu_int)^2,
                )),
                energy = energy_bc,
                moisture = moisture_bc,
                turbconv = turbconv_bcs(turbconv),
            ),
            AtmosBC(),
        ),
        init_state_prognostic = ics,
    )

    # Assemble model components
    model = AtmosModel{FT}(
        config_type,
        param_set;
        problem = problem,
        turbulence = SmagorinskyLilly{FT}(C_smag),
        moisture = moisture,
        source = source,
        turbconv = turbconv,
    )

    return model
end

function config_diagnostics(driver_config)
    default_dgngrp = setup_atmos_default_diagnostics(
        AtmosLESConfigType(),
        "2500steps",
        driver_config.name,
    )
    core_dgngrp = setup_atmos_core_diagnostics(
        AtmosLESConfigType(),
        "2500steps",
        driver_config.name,
    )
    return ClimateMachine.DiagnosticsConfiguration([
        default_dgngrp,
        core_dgngrp,
    ])
end
