#!/usr/bin/env julia --project
#=
# This experiment file establishes the initial conditions, boundary conditions,
# source terms and simulation parameters (domain size + resolution) for the
# gabls LES case. The set of parameters presented in the `master` branch copy
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
# julia --project experiments/AtmosLES/gabls.jl
#
# See `src/Driver/driver_configs.jl` for additional flags (e.g. VTK, diagnostics,
# update-interval, output directory settings)
#
# Upcoming changes:
# 1) Atomic sources
# 2) Improved boundary conditions
# 3) Collapsed experiment design
# 4) Updates to generally keep this in sync with master

Beare, R.J., Macvean, M.K., Holtslag, A.A.M. et al.
An Intercomparison of Large-Eddy Simulations of the Stable Boundary Layer.
Boundary-Layer Meteorol 118, 247–272 (2006). https://doi.org/10.1007/s10546-004-2820-6
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
  gabls Geostrophic Forcing (Source)
"""
struct GablsGeostrophic{FT} <: Source
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
    s::GablsGeostrophic,
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
  gabls Sponge (Source)
"""
struct GablsSponge{FT} <: Source
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
    s::GablsSponge,
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
  GablsTendencies (Source)
No tendencies are applied in gabls
"""

"""
  Initial Condition for gabls LES
"""
function init_gabls!(problem, bl, state, aux, (x, y, z), t)
    # This experiment runs the gabls LES Configuration
    # (Stable boundary layer regime)
    # x,y,z imply eastward, northward and altitude coordinates in `[m]`

    # Problem floating point precision
    FT = eltype(state)

    P_sfc::FT = 1.0e5 # Surface air pressuref
    _R_d::FT = R_d(bl.param_set)
    θ_dry_sfc = FT(265) # Prescribed θ_dry at surface
    T_sfc = FT(265) # Surface temperature
    _grav = FT(grav(bl.param_set))

    # Initialise speeds [u = Eastward, v = Northward, w = Vertical]
    u::FT = 8.0
    v::FT = 0
    w::FT = 0

    # Assign piecewise quantities to θ_dry
    θ_dry::FT = 0

    # Piecewise functions for potential temperature and total moisture
    if FT(0) <= z <= FT(100)
        # Well mixed layer
        θ_dry = FT(265)
    else
        # Stable layer
        θ_dry = FT(265) + (z - FT(100.0)) * FT(0.01)
    end
    # Scale height based on surface parameters
    H = _R_d * T_sfc / _grav
    # Pressure based on scale height
    P = P_sfc * exp(-z / H)

    # Establish thermodynamic state and dry phase partitioning
    TS = PhaseDry_given_pθ(bl.param_set, P, θ_dry)
    T = air_temperature(TS)
    ρ = air_density(TS)

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

    if z <= FT(100) # Add random perturbations to bottom 100m of model
        state.ρe += (rand() - FT(0.5)) * ρe_tot / 100
    end
    init_state_prognostic!(bl.turbconv, bl, state, aux, (x, y, z), t)
end

function gabls_model(
    ::Type{FT},
    config_type,
    zmax,
    surface_flux;
    turbconv = NoTurbConv(),
) where {FT}

    ics = init_gabls!     # Initial conditions

    C_smag = FT(0.23)     # Smagorinsky coefficient

    u_star = FT(0.25)     # Friction velocity
    C_drag = FT(0.0011)   # Bulk transfer coefficient

    T_sfc = FT(265)     # Surface temperature `[K]`
    SHF = FT(0)         # Sensible heat flux `[W/m²]`

    z_sponge = FT(100)     # Start of sponge layer
    α_max = FT(0.75)        # Strength of sponge layer (timescale)
    γ = 2              # Strength of sponge layer (exponent)

    u_geostrophic = FT(8)        # Eastward relaxation speed
    u_slope = FT(0)     # Slope of altitude-dependent relaxation speed
    v_geostrophic = FT(0)          # Northward relaxation speed
    f_coriolis = FT(1.39e-4) # Coriolis parameter (latitude 73 N)

    # Assemble source components
    source_default = (
        Gravity(),
        GablsSponge{FT}(
            zmax,
            z_sponge,
            α_max,
            γ,
            u_geostrophic,
            u_slope,
            v_geostrophic,
        ),
        GablsGeostrophic{FT}(f_coriolis, u_geostrophic, u_slope, v_geostrophic),
        turbconv_sources(turbconv)...,
    )

    source = source_default
    # Set up problem initial and boundary conditions
    if surface_flux == "prescribed"
        energy_bc = PrescribedTemperature((state, aux, t) -> T_sfc - FT(0.25/3600.0)*t)
    elseif surface_flux == "bulk"
        energy_bc = BulkFormulaEnergy(
            (state, aux, t, normPu_int) -> C_drag,
            (state, aux, t) -> (T_sfc),
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
        moisture = DryModel(),
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
