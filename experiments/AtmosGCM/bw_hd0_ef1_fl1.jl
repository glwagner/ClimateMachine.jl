#!/usr/bin/env julia --project
using ClimateMachine
ClimateMachine.cli()

using ClimateMachine.Atmos
using ClimateMachine.ConfigTypes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.Orientations
using ClimateMachine.ODESolvers
using ClimateMachine.SystemSolvers: ManyColumnLU
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.TemperatureProfiles
using ClimateMachine.Thermodynamics: air_density, total_energy
using ClimateMachine.VariableTemplates
using ClimateMachine.TurbulenceClosures

using LinearAlgebra
using StaticArrays
using Test

using CLIMAParameters
using CLIMAParameters.Planet: MSLP, R_d, day, grav, Omega, planet_radius
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

function init_baroclinic_wave!(bl, state, aux, coords, t)
    FT = eltype(state)

    # parameters 
    _grav::FT = grav(bl.param_set)
    _R_d::FT = R_d(bl.param_set)
    _Ω::FT = Omega(bl.param_set)
    _a::FT = planet_radius(bl.param_set)
    _p_0::FT = MSLP(bl.param_set)
    
    k::FT = 3
    T_E::FT = 310
    T_P::FT = 240
    T_0::FT = 0.5 * (T_E + T_P)
    Γ::FT = 0.005
    A::FT = 1 / Γ
    B::FT = (T_0-T_P)/T_0/T_P
    C::FT = 0.5 * (k+2) * (T_E-T_P)/T_E/T_P
    b::FT = 2
    H::FT = _R_d * T_0 / _grav
    z_t::FT = 15e3
    λ_c::FT = π/9
    φ_c::FT = 2*π/9
    d_0::FT = _a/6
    V_p::FT = 1

    # grid
    φ = latitude(bl, aux)
    λ = longitude(bl, aux)
    z = altitude(bl, aux)
    r::FT = z+_a
    γ::FT = 0 # set to 0 for shallow-atmosphere case and to 1 for deep atmosphere case

    # convenience functions for temperature and pressure
    τ_z_1::FT = exp(Γ*z/T_0)
    τ_z_2::FT = 1 - 2*(z/b/H)^2
    τ_z_3::FT = exp(-(z/b/H)^2)
    τ_1::FT = 1/T_0 * τ_z_1 + B * τ_z_2 * τ_z_3 
    τ_2::FT = C * τ_z_2 * τ_z_3 
    τ_int_1::FT = A * (τ_z_1-1) + B * z * τ_z_3
    τ_int_2::FT = C * z * τ_z_3
    I_T::FT = (cos(φ) * (1 + γ*z/_a))^k - k/(k+2) * (cos(φ) * (1 + γ*z/_a))^(k+2) 

    # base state temperature, pressure, specific humidity, density
    T::FT = (1/(1 + γ*z/_a))^2 * (τ_1 - τ_2 * I_T)^(-1) 
    p::FT = _p_0 * exp(-_grav/_R_d * (τ_int_1 - τ_int_2 * I_T))
    
    # base state velocity
    U::FT = _grav*k/_a * τ_int_2 * T * ((cos(φ) * (1 + γ*z/_a))^(k-1) - (cos(φ) * (1 + γ*z/_a))^(k+1))
    u_ref::FT = -_Ω*(_a + γ*z) *cos(φ) + sqrt((_Ω*(_a + γ*z)*cos(φ))^2 + (_a + γ*z)*cos(φ)*U)
    v_ref::FT = 0
    w_ref::FT = 0

    # velocity perturbations
    F_z::FT = 1 - 3*(z/z_t)^2 + 2*(z/z_t)^3 
    if z > z_t
      F_z = FT(0)
    end
    d::FT = _a*acos(sin(φ)*sin(φ_c) + cos(φ)*cos(φ_c)*cos(λ-λ_c)) 
    c3::FT = cos(π*d/2/d_0)^3
    s1::FT = sin(π*d/2/d_0)
    if 0 < d < d_0 && d != FT(_a*π)
      u′::FT = -16*V_p/3/sqrt(3) * F_z * c3 * s1 * (-sin(φ_c)*cos(φ) + cos(φ_c)*sin(φ)*cos(λ-λ_c)) / sin(d/_a) 
      v′::FT = 16*V_p/3/sqrt(3) * F_z * c3 * s1 * cos(φ_c)*sin(λ-λ_c) / sin(d/_a)
    else
      u′ = FT(0)
      v′ = FT(0)
    end
    w′::FT = 0
    u_sphere = SVector{3, FT}(u_ref + u′, v_ref + v′, w_ref + w′)
    u_cart = sphr_to_cart_vec(bl.orientation, u_sphere, aux)
    
    # potential & kinetic energy
    ρ = air_density(bl.param_set, T, p)
    e_pot::FT = gravitational_potential(bl.orientation, aux)
    e_kin::FT = 0.5 * u_cart' * u_cart
    e_tot::FT = total_energy(bl.param_set, e_kin, e_pot, T)

    state.ρ = ρ
    state.ρu = ρ * u_cart
    state.ρe = ρ * e_tot
    
    nothing
end

function config_baroclinic_wave(FT, poly_order, resolution)
    # Set up a reference state for linearization of equations
    temp_profile_ref = DecayingTemperatureProfile{FT}(param_set, FT(275), FT(75), FT(45e3))
    ref_state = HydrostaticState(temp_profile_ref)
    domain_height::FT = 30e3 # distance between surface and top of atmosphere (m)
    
    # Set up the atmosphere model
    exp_name = "BaroclinicWave"

    model = AtmosModel{FT}(
        AtmosGCMConfigType,
        param_set;
        ref_state = ref_state,
        turbulence = ConstantViscosityWithDivergence(FT(0.0)),
        hyperdiffusion = NoHyperDiffusion(),
        #hyperdiffusion = StandardHyperDiffusion(FT(1000)),
        moisture = DryModel(),
        source = (Gravity(), Coriolis(),),
        init_state_conservative = init_baroclinic_wave!,
    )

    config = ClimateMachine.AtmosGCMConfiguration(
        exp_name,
        poly_order,
        resolution,
        domain_height,
        param_set,
        init_baroclinic_wave!;
        model = model,
    )

    return config
end

function main()
    # Driver configuration parameters
    FT = Float64                             # floating type precision
    poly_order = 4                           # discontinuous Galerkin polynomial order
    n_horz = 7                               # horizontal element number
    n_vert = 2                               # vertical element number
    n_days::FT = 15
    timestart::FT = 0                        # start time (s)
    timeend::FT = n_days * day(param_set)    # end time (s)

    # Set up driver configuration
    driver_config = config_baroclinic_wave(FT, poly_order, (n_horz, n_vert))

    ode_solver_type = ClimateMachine.ExplicitSolverType(
         solver_method = LSRK144NiegemannDiehlBusch,
    )
   # ode_solver_type = ClimateMachine.IMEXSolverType(
   #     implicit_model = AtmosAcousticGravityLinearModel,
   #     implicit_solver = ManyColumnLU,
   #     solver_method = ARK2GiraldoKellyConstantinescu,
   #     split_explicit_implicit = true,
   #     discrete_splitting = false,
   # )
    
    # Set up experiment
    CFL::FT = 0.1
    solver_config = ClimateMachine.SolverConfiguration(
        timestart,
        timeend,
        driver_config,
        Courant_number = CFL,
        ode_solver_type = ode_solver_type,
#        CFL_direction = HorizontalDirection(),
        CFL_direction = EveryDirection(),
#        diffdir= HorizontalDirection(),
    )

    # Set up diagnostics
    dgn_config = config_diagnostics(FT, driver_config)

    # Set up user-defined callbacks
    filterorder = 32
    filter = ExponentialFilter(solver_config.dg.grid, 0, filterorder)
    cbfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            AtmosFilterPerturbations(driver_config.bl),
            solver_config.dg.grid,
            filter,
            state_auxiliary = solver_config.dg.state_auxiliary,
        )
        nothing
    end

    # Run the model
    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (cbfilter,),
        check_euclidean_distance = true,
    )
end

function config_diagnostics(FT, driver_config)
    interval = "40000steps" # chosen to allow a single diagnostics collection

    _planet_radius = FT(planet_radius(param_set))

    info = driver_config.config_info
    boundaries = [
        FT(-90.0) FT(-180.0) _planet_radius
        FT(90.0) FT(180.0) FT(_planet_radius + info.domain_height)
    ]
    resolution = (FT(1), FT(1), FT(500)) # in (deg, deg, m)
    interpol = ClimateMachine.InterpolationConfiguration(
        driver_config,
        boundaries,
        resolution,
    )

    dgngrp = setup_atmos_default_diagnostics(
        AtmosGCMConfigType(),
        interval,
        driver_config.name,
        interpol = interpol,
    )

    return ClimateMachine.DiagnosticsConfiguration([dgngrp])
end

main()
