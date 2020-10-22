using ClimateMachine
ClimateMachine.init(;
    parse_clargs = true,
    output_dir = get(ENV, "CLIMATEMACHINE_SETTINGS_OUTPUT_DIR", "output"),
    fix_rng_seed = true,
)
using ClimateMachine.SingleStackUtils
using ClimateMachine.Checkpoint
using ClimateMachine.BalanceLaws: vars_state
const clima_dir = dirname(dirname(pathof(ClimateMachine)));

include(joinpath(clima_dir, "experiments", "AtmosLES", "bomex_model.jl"))
include("edmf_model.jl")
include("edmf_kernels.jl")

"""
    init_state_prognostic!(
            turbconv::EDMF{FT},
            m::AtmosModel{FT},
            state::Vars,
            aux::Vars,
            coords,
            t::Real,
        ) where {FT}

Initialize EDMF state variables.
This method is only called at `t=0`.
"""
function init_state_prognostic!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
) where {FT}
    # Aliases:
    gm = state
    en = state.turbconv.environment
    up = state.turbconv.updraft
    N_up = n_updrafts(turbconv)
    # GCM setting - Initialize the grid mean profiles of prognostic variables (ρ,e_int,q_tot,u,v,w)
    z = altitude(m, aux)

    # SCM setting - need to have separate cases coded and called from a folder - see what LES does
    # a moist_thermo state is used here to convert the input θ,q_tot to e_int, q_tot profile
    e_int = internal_energy(m, state, aux)

    q_tot = state.moisture.ρq_tot / state.ρ
    ts = PhaseEquil(m.param_set, e_int, state.ρ, q_tot)
    T = air_temperature(ts)
    p = air_pressure(ts)
    q = PhasePartition(ts)
    θ_liq = liquid_ice_pottemp(ts)

    a_min = turbconv.subdomains.a_min
    @unroll_map(N_up) do i
        up[i].ρa = gm.ρ * a_min
        up[i].ρaw = gm.ρu[3] * a_min
        up[i].ρaθ_liq = gm.ρ * a_min * θ_liq
        up[i].ρaq_tot = gm.moisture.ρq_tot * a_min
    end

    # initialize environment covariance with zero for now
    if z <= FT(2500)
        en.ρatke = gm.ρ * (FT(1) - z / FT(3000))
    else
        en.ρatke = FT(0)
    end
    en.ρaθ_liq_cv = FT(1e-5) / max(z, FT(10))
    en.ρaq_tot_cv = FT(1e-5) / max(z, FT(10))
    en.ρaθ_liq_q_tot_cv = FT(1e-7) / max(z, FT(10))
    return nothing
end;

function main(::Type{FT}) where {FT}
    # add a command line argument to specify the kind of surface flux
    # TODO: this will move to the future namelist functionality
    bomex_args = ArgParseSettings(autofix_names = true)
    add_arg_group!(bomex_args, "BOMEX")
    @add_arg_table! bomex_args begin
        "--surface-flux"
        help = "specify surface flux for energy and moisture"
        metavar = "prescribed|bulk"
        arg_type = String
        default = "prescribed"
    end

    cl_args =
        ClimateMachine.init(parse_clargs = true, custom_clargs = bomex_args)

    surface_flux = cl_args["surface_flux"]

    # DG polynomial order
    N = 1
    nelem_vert = 50

    # Prescribe domain parameters
    zmax = FT(3000)

    t0 = FT(0)

    # Simulation time
    timeend = FT(6*3600)
    # timeend = FT(1000)
    CFLmax = FT(0.90)

    config_type = SingleStackConfigType

    ode_solver_type = ClimateMachine.ExplicitSolverType(
        solver_method = LSRK144NiegemannDiehlBusch,
    )

    N_updrafts = 1
    N_quad = 3
    turbconv = EDMF(FT, N_updrafts, N_quad, zmax)

    model =
        bomex_model(FT, config_type, zmax, surface_flux; turbconv = turbconv)

    # Assemble configuration
    driver_config = ClimateMachine.SingleStackConfiguration(
        "BOMEX_EDMF",
        N,
        nelem_vert,
        zmax,
        param_set,
        model;
        hmax = FT(500),
        solver_type = ode_solver_type,
    )

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        init_on_cpu = true,
        Courant_number = CFLmax,
        # fixed_number_of_steps = 1000,
        # fixed_number_of_steps=1082 # last timestep before crash
    )

    # --- Zero-out horizontal variations:
    vsp = vars_state(model, Prognostic(), FT)
    horizontally_average!(
        driver_config.grid,
        solver_config.Q,
        varsindex(vsp, :turbconv),
    )
    horizontally_average!(
        driver_config.grid,
        solver_config.Q,
        varsindex(vsp, :ρe),
    )
    horizontally_average!(
        driver_config.grid,
        solver_config.Q,
        varsindex(vsp, :moisture, :ρq_tot),
    )

    vsa = vars_state(model, Auxiliary(), FT)
    horizontally_average!(
        driver_config.grid,
        solver_config.dg.state_auxiliary,
        varsindex(vsa, :turbconv),
    )
    # ---

    dgn_config = config_diagnostics(driver_config)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            ("moisture.ρq_tot", turbconv_filters(turbconv)...),
            solver_config.dg.grid,
            TMARFilter(),
        )
        bl = solver_config.dg.balance_law
        st = vars_state(bl, Prognostic(), eltype(Q))
        aux = vars_state(bl, Auxiliary(), eltype(Q))
        a_min = bl.turbconv.subdomains.a_min
        a_max = bl.turbconv.subdomains.a_max
        i_ρ = varsindex(st, :ρ)
        i_ρq_tot_gm = varsindex(st, :moisture, :ρq_tot)
        ρ_gm = solver_config.Q[:,i_ρ,:]
        i_ρaθ_liq_en = varsindex(aux, :turbconv, :environment,:ρaθ_liq)
        ρaθ_liq_en  = solver_config.dg.state_auxiliary[:,i_ρaθ_liq_en,:]
        ρaθ_liq_ups = FT(0)
        ρaq_tot_ups = FT(0)
        ρa_ups      = FT(0)
        ρaw_ups     = FT(0)
        for i in 1:n_updrafts(bl.turbconv)
            i_ρaθ_liq_up = varsindex(st, :turbconv, :updraft, Val(i), :ρaθ_liq)
            i_ρaq_tot_up = varsindex(st, :turbconv, :updraft, Val(i), :ρaq_tot)
            i_ρa_up      = varsindex(st, :turbconv, :updraft, Val(i), :ρa)
            i_ρaw_up     = varsindex(st, :turbconv, :updraft, Val(i), :ρaw)
            ρaθ_liq_ups  = ρaθ_liq_ups .+ solver_config.Q[:,i_ρaθ_liq_up,:]
            ρaq_tot_ups  = ρaq_tot_ups .+ solver_config.Q[:,i_ρaq_tot_up,:]
            ρa_ups       = ρa_ups .+ solver_config.Q[:,i_ρa_up,:]
            ρaw_ups      = ρaw_ups .+ solver_config.Q[:,i_ρaw_up,:]
        end
        ρq_tot_gm   = solver_config.Q[:,i_ρq_tot_gm,:]

        ρaw_en      = - ρaw_ups
        ρaq_tot_en  = (ρq_tot_gm - ρaq_tot_ups) ./ (ρ_gm .- ρa_ups)
        θ_liq_en    = ρaθ_liq_en ./ (ρ_gm .- ρa_ups)
        q_tot_en    = ρaq_tot_en ./ (ρ_gm .- ρa_ups)
        w_en        = ρaw_en./ (ρ_gm .- ρa_ups)
        for i in 1:n_updrafts(bl.turbconv)
            i_ρa_up      = varsindex(st, :turbconv, :updraft, Val(i), :ρa)
            i_ρaθ_liq_up = varsindex(st, :turbconv, :updraft, Val(i), :ρaθ_liq)
            i_ρaq_tot_up = varsindex(st, :turbconv, :updraft, Val(i), :ρaq_tot)
            i_ρaw_up     = varsindex(st, :turbconv, :updraft, Val(i), :ρaw)

            ρa_up_mask_min   = solver_config.Q[:,i_ρa_up,:] .< ρ_gm*a_min
            ρa_up_mask_max   = solver_config.Q[:,i_ρa_up,:] .> ρ_gm*a_max
            ρa_up    = solver_config.Q[:,i_ρa_up,:]
            ρ_δa_min = max.(ρa_up_mask_min.*(ρ_gm*a_min.-ρa_up),FT(0))
            ρ_δa_max = min.(ρa_up_mask_max.*(ρ_gm*a_max.-ρa_up),FT(0))

            solver_config.Q[:,i_ρa_up,:]      .+= ρa_up_mask_min.*ρ_δa_min
            solver_config.Q[:,i_ρaθ_liq_up,:] .+= ρa_up_mask_min.*ρ_δa_min.*θ_liq_en
            solver_config.Q[:,i_ρaq_tot_up,:] .+= ρa_up_mask_min.*ρ_δa_min.*q_tot_en
            solver_config.Q[:,i_ρaw_up,:]     .+= ρa_up_mask_min.*ρ_δa_min.*w_en

            ρa_up    = solver_config.Q[:,i_ρa_up,:]
            θ_liq_up = solver_config.Q[:,i_ρaθ_liq_up,:] ./ ρa_up
            q_tot_up = solver_config.Q[:,i_ρaq_tot_up,:] ./ ρa_up
            w_up     = solver_config.Q[:,i_ρaw_up,:]     ./ ρa_up

            solver_config.Q[:,i_ρa_up,:]      .+= ρa_up_mask_max.*ρ_δa_max
            solver_config.Q[:,i_ρaθ_liq_up,:] .+= ρa_up_mask_max.*ρ_δa_max.*θ_liq_up
            solver_config.Q[:,i_ρaq_tot_up,:] .+= ρa_up_mask_max.*ρ_δa_max.*q_tot_up
            solver_config.Q[:,i_ρaw_up,:]     .+= ρa_up_mask_max.*ρ_δa_max.*w_up

        end
        nothing
    end

    # State variable
    Q = solver_config.Q
    # Volume geometry information
    vgeo = driver_config.grid.vgeo
    M = vgeo[:, Grids._M, :]
    # Unpack prognostic vars
    ρ₀ = Q.ρ
    ρe₀ = Q.ρe
    # DG variable sums
    Σρ₀ = sum(ρ₀ .* M)
    Σρe₀ = sum(ρe₀ .* M)

    grid = driver_config.grid

    # state_types = (Prognostic(), Auxiliary(), GradientFlux())
    state_types = (Prognostic(), Auxiliary())
    all_data = [dict_of_nodal_states(solver_config, ["z"], state_types)]
    time_data = FT[0]

    # Define the number of outputs from `t0` to `timeend`
    n_outputs = 10
    # This equates to exports every ceil(Int, timeend/n_outputs) time-step:
    every_x_simulation_time = ceil(Int, timeend / n_outputs)

    cb_data_vs_time =
        GenericCallbacks.EveryXSimulationTime(every_x_simulation_time) do
            push!(
                all_data,
                dict_of_nodal_states(solver_config, ["z"], state_types),
            )
            push!(time_data, gettime(solver_config.solver))
            nothing
        end

    cb_check_cons = GenericCallbacks.EveryXSimulationSteps(3000) do
        Q = solver_config.Q
        δρ = (sum(Q.ρ .* M) - Σρ₀) / Σρ₀
        δρe = (sum(Q.ρe .* M) .- Σρe₀) ./ Σρe₀
        @show (abs(δρ))
        @show (abs(δρe))
        @test (abs(δρ) <= 0.001)
        @test (abs(δρe) <= 0.0025)
        nothing
    end

    cb_print_step = GenericCallbacks.EveryXSimulationSteps(100) do
        @show getsteps(solver_config.solver)
        nothing
    end

    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (
            cbtmarfilter,
            cb_check_cons,
            cb_data_vs_time,
            cb_print_step,
        ),
        check_euclidean_distance = true,
    )

    push!(all_data, dict_of_nodal_states(solver_config, ["z"], state_types))
    push!(time_data, gettime(solver_config.solver))

    return solver_config, all_data, time_data, state_types
end

solver_config, all_data, time_data, state_types = main(Float64)

include(joinpath(
    clima_dir,
    "test",
    "Atmos",
    "EDMF",
    "bomex_edmf_regression_test.jl",
))
