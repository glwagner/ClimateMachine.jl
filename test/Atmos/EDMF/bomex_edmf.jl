using ClimateMachine
using ClimateMachine.SingleStackUtils
using ClimateMachine.BalanceLaws: vars_state
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
using Test
include(joinpath(clima_dir, "experiments", "AtmosLES", "bomex_model.jl"))
include("edmf_model.jl")
include("edmf_kernels.jl")

# ===============================================
# case configuration / flow control
struct UseExplicitSolver end
struct UseImplicitSolver end

struct InitWithZero end
struct InitWithNaN end

get_tc_data(::InitWithZero) = 0
get_tc_data(::InitWithNaN) = NaN
const cases = (
    (solver_type=UseExplicitSolver(),turbconv_data=InitWithZero()), # works
    (solver_type=UseExplicitSolver(),turbconv_data=InitWithNaN()),  # works
    (solver_type=UseImplicitSolver(),turbconv_data=InitWithZero()), # works (with ode_dt prescribed)
    (solver_type=UseImplicitSolver(),turbconv_data=InitWithNaN()),  # fails
    )
# ===============================================

function init_state_prognostic!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
) where {FT}

    N_up = n_updrafts(turbconv)
    val = get_tc_data(turbconv.turbconv_data)
    for i in 1:N_up
        state.turbconv.updraft[i].ρa = val
    end
    return nothing
end;
# ===============================================

function main(::Type{FT}, solver_type, turbconv_data) where {FT}
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
    timeend = FT(100)
    CFLmax = FT(0.90)
    config_type = SingleStackConfigType

    if solver_type isa UseImplicitSolver
        ode_solver_type = ClimateMachine.IMEXSolverType()
    elseif solver_type isa UseExplicitSolver
        ode_solver_type = ClimateMachine.ExplicitSolverType(
            solver_method = LSRK144NiegemannDiehlBusch,
        )
    else
        error("Bad option")
    end

    N_updrafts = 1
    turbconv = EDMF(FT, N_updrafts; turbconv_data=turbconv_data)

    model = bomex_model(FT, config_type, zmax, surface_flux, turbconv)

    # Assemble configuration
    driver_config = ClimateMachine.SingleStackConfiguration(
        "BOMEX_EDMF",
        N,
        nelem_vert,
        zmax,
        param_set,
        model;
        solver_type = ode_solver_type,
    )

    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        init_on_cpu = true,
        ode_dt = FT(1.54979e-01),
        Courant_number = CFLmax,
    )

    dgn_config = config_diagnostics(driver_config)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            (
                "moisture.ρq_tot",
                turbconv_filters(turbconv)...,
            ),
            solver_config.dg.grid,
            TMARFilter(),
        )
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

    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (cbtmarfilter, cb_check_cons),
        check_euclidean_distance = true,
    )

    all_data = [dict_of_nodal_states(solver_config, ["z"], (Prognostic(), Auxiliary()))]
    return all_data
end

# We can comment out `ode_dt = FT(1.54979e-01)` above for some combination
@testset "Testing Explicit/Implicit solvers" begin
    for i in 1:4
        solver_type = cases[i].solver_type # used above
        turbconv_data = cases[i].turbconv_data # used above
        @show solver_type, turbconv_data

        all_data = main(Float64, solver_type, turbconv_data)
        @test !isnan(norm(all_data[end]["ρ"]))
        @test !isnan(norm(all_data[end]["ρu[1]"]))
        @test !isnan(norm(all_data[end]["ρu[2]"]))
        @test !isnan(norm(all_data[end]["ρu[3]"]))
        @test !isnan(norm(all_data[end]["ρe"]))
    end
end
