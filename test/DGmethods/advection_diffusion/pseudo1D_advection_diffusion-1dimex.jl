using MPI
using CLIMA
using Logging
using Test
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.DGmethods
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.LinearSolvers
using CLIMA.GeneralizedMinimalResidualSolver
using CLIMA.ColumnwiseLUSolver: SingleColumnLU, ManyColumnLU, banded_matrix,
                                banded_matrix_vector_product!
using CLIMA.AdditiveRungeKuttaMethod
using LinearAlgebra
using Printf
using Dates
using CLIMA.GenericCallbacks: EveryXWallTimeSeconds, EveryXSimulationSteps
using CLIMA.ODESolvers: solve!, gettime
using CLIMA.VTK: writevtk, writepvtu
using CLIMA.DGmethods: EveryDirection, HorizontalDirection, VerticalDirection

@static if haspkg("CuArrays")
  using CUDAdrv
  using CUDAnative
  using CuArrays
  CuArrays.allowscalar(false)
  const ArrayTypes = (CuArray, )
else
  const ArrayTypes = (Array, )
end

if !@isdefined integration_testing
  if length(ARGS) > 0
    const integration_testing = parse(Bool, ARGS[1])
  else
    const integration_testing =
      parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
  end
end

const output = parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_OUTPUT","false")))

include("advection_diffusion_model.jl")

struct Pseudo1D{n, α, β, μ, δ} <: AdvectionDiffusionProblem end

function init_velocity_diffusion!(::Pseudo1D{n, α, β}, aux::Vars,
                                  geom::LocalGeometry) where {n, α, β}
  # Direction of flow is n with magnitude α
  aux.u = α * n

  # diffusion of strength β in the n direction
  aux.D = β * n * n'
end

function initial_condition!(::Pseudo1D{n, α, β, μ, δ}, state, aux, x,
                            t) where {n, α, β, μ, δ}
  ξn = dot(n, x)
  # ξT = SVector(x) - ξn * n
  state.ρ = exp(-(ξn - μ - α * t)^2 / (4 * β * (δ + t))) / sqrt(1 + t / δ)
end

function do_output(mpicomm, vtkdir, vtkstep, dg, Q, Qe, model, testname)
  ## name of the file that this MPI rank will write
  filename = @sprintf("%s/%s_mpirank%04d_step%04d",
                      vtkdir, testname, MPI.Comm_rank(mpicomm), vtkstep)

  statenames = flattenednames(vars_state(model, eltype(Q)))
  exactnames = statenames .* "_exact"

  writevtk(filename, Q, dg, statenames, Qe, exactnames)

  ## Generate the pvtu file for these vtk files
  if MPI.Comm_rank(mpicomm) == 0
    ## name of the pvtu file
    pvtuprefix = @sprintf("%s/%s_step%04d", vtkdir, testname, vtkstep)

    ## name of each of the ranks vtk files
    prefixes = ntuple(MPI.Comm_size(mpicomm)) do i
      @sprintf("%s_mpirank%04d_step%04d", testname, i - 1, vtkstep)
    end

    writepvtu(pvtuprefix, prefixes, (statenames..., exactnames...))

    @info "Done writing VTK: $pvtuprefix"
  end
end


function run(mpicomm, ArrayType, dim, topl, N, timeend, FT, dt,
             n, α, β, μ, δ, vtkdir, outputtime, linearsolvertype)

  grid = DiscontinuousSpectralElementGrid(topl,
                                          FloatType = FT,
                                          DeviceArray = ArrayType,
                                          polynomialorder = N,
                                         )
  model = AdvectionDiffusion{dim}(Pseudo1D{n, α, β, μ, δ}())
  dg = DGModel(model,
               grid,
               Rusanov(),
               CentralNumericalFluxDiffusive(),
               CentralGradPenalty(),
               direction=EveryDirection())

  hdg = DGModel(model,
                grid,
                Rusanov(),
                CentralNumericalFluxDiffusive(),
                CentralGradPenalty(),
                auxstate=dg.auxstate,
                direction=HorizontalDirection())

  vdg = DGModel(model,
                grid,
                Rusanov(),
                CentralNumericalFluxDiffusive(),
                CentralGradPenalty(),
                auxstate=dg.auxstate,
                direction=VerticalDirection())


  let
    # test that we can build the banded matrix
    A_banded =
      banded_matrix(vdg, MPIStateArray(dg), MPIStateArray(dg);
                    single_column=linearsolvertype <: SingleColumnLU)

    Q = MPIStateArray(vdg)
    dQ1 = MPIStateArray(vdg)
    dQ2 = MPIStateArray(vdg)
    Q.data .= rand(size(Q.data))
    vdg(dQ1, Q, nothing, 0; increment=false)

    banded_matrix_vector_product!(vdg, A_banded, dQ2, Q)
    @test dQ1.realdata ≈ dQ2.realdata
  end

  Q = init_ode_state(dg, FT(0))

  ode_solver = ARK548L2SA2KennedyCarpenter(hdg, vdg, linearsolvertype(), Q;
                                           dt = dt, t0 = 0,
                                           split_nonlinear_linear=true)


  eng0 = norm(Q)
  @info @sprintf """Starting
  norm(Q₀) = %.16e""" eng0

  # Set up the information callback
  starttime = Ref(now())
  cbinfo = EveryXWallTimeSeconds(60, mpicomm) do (s=false)
    if s
      starttime[] = now()
    else
      energy = norm(Q)
      @info @sprintf("""Update
                     simtime = %.16e
                     runtime = %s
                     norm(Q) = %.16e""", gettime(ode_solver),
                     Dates.format(convert(Dates.DateTime,
                                          Dates.now()-starttime[]),
                                  Dates.dateformat"HH:MM:SS"),
                     energy)
    end
  end
  callbacks = (cbinfo,)
  if ~isnothing(vtkdir)
    # create vtk dir
    mkpath(vtkdir)

    vtkstep = 0
    # output initial step
    do_output(mpicomm, vtkdir, vtkstep, dg, Q, Q, model, "advection_diffusion")

    # setup the output callback
    cbvtk = EveryXSimulationSteps(floor(outputtime/dt)) do
      vtkstep += 1
      Qe = init_ode_state(dg, gettime(ode_solver))
      do_output(mpicomm, vtkdir, vtkstep, dg, Q, Qe, model,
                "advection_diffusion")
    end
    callbacks = (callbacks..., cbvtk)
  end

  numberofsteps = convert(Int64, cld(timeend, dt))
  dt = timeend / numberofsteps

  @info "time step" dt numberofsteps dt*numberofsteps timeend

  solve!(Q, ode_solver; numberofsteps=numberofsteps, callbacks=callbacks,
         adjustfinalstep=false)

  # Print some end of the simulation information
  engf = norm(Q)
  Qe = init_ode_state(dg, FT(timeend))

  engfe = norm(Qe)
  errf = euclidean_distance(Q, Qe)
  @info @sprintf """Finished
  norm(Q)                 = %.16e
  norm(Q) / norm(Q₀)      = %.16e
  norm(Q) - norm(Q₀)      = %.16e
  norm(Q - Qe)            = %.16e
  norm(Q - Qe) / norm(Qe) = %.16e
  """ engf engf/eng0 engf-eng0 errf errf / engfe
  errf
end

let
  MPI.Initialized() || MPI.Init()
  mpicomm = MPI.COMM_WORLD
  ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
  loglevel = ll == "DEBUG" ? Logging.Debug :
  ll == "WARN"  ? Logging.Warn  :
  ll == "ERROR" ? Logging.Error : Logging.Info
  logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
  global_logger(ConsoleLogger(logger_stream, loglevel))
  @static if haspkg("CUDAnative")
    device!(MPI.Comm_rank(mpicomm) % length(devices()))
  end

  polynomialorder = 4
  base_num_elem = 4

  expected_result = Dict()
  expected_result[2, 1, Float64] = 1.2228434091602128e-02
  expected_result[2, 2, Float64] = 8.8037798002260420e-04
  expected_result[2, 3, Float64] = 4.8828676920661276e-05
  expected_result[2, 4, Float64] = 2.0105646643725454e-06
  expected_result[3, 1, Float64] = 9.5425450102548364e-03
  expected_result[3, 2, Float64] = 5.9769045240778518e-04
  expected_result[3, 3, Float64] = 4.0081798525590592e-05
  expected_result[3, 4, Float64] = 2.9803558844543670e-06
  expected_result[2, 1, Float32] = 1.2228445149958134e-02
  expected_result[2, 2, Float32] = 8.8042858988046646e-04
  expected_result[2, 3, Float32] = 4.8848683945834637e-05
  expected_result[2, 4, Float32] = 2.1847356492799008e-06
  expected_result[3, 1, Float32] = 9.5424978062510490e-03
  expected_result[3, 2, Float32] = 5.9770536608994007e-04
  expected_result[3, 3, Float32] = 4.0205955883720890e-05
  expected_result[3, 4, Float32] = 5.1591650844784454e-06

  numlevels = integration_testing ? 4 : 1

  @testset "$(@__FILE__)" for ArrayType in ArrayTypes
    for FT in (Float64, Float32)
      result = zeros(FT, numlevels)
      for dim = 2:3
        for linearsolvertype in (SingleColumnLU, ManyColumnLU)
          d = dim == 2 ? FT[1, 10, 0] : FT[1, 1, 10]
          n = SVector{3, FT}(d ./ norm(d))

          α = FT(1)
          β = FT(1 // 100)
          μ = FT(-1 // 2)
          δ = FT(1 // 10)
          for l = 1:numlevels
            Ne = 2^(l-1) * base_num_elem
            brickrange = (ntuple(j->range(FT(-1); length=Ne+1, stop=1), dim-1)...,
                          range(FT(-5); length=5Ne+1, stop=5))

            periodicity = ntuple(j->false, dim)
            topl = StackedBrickTopology(mpicomm, brickrange;
                                        periodicity = periodicity,
                                        boundary = (ntuple(j->(1,1), dim-1)...,
                                                    (3,3)))
            dt = (α/4) / (Ne * polynomialorder^2)

            outputtime = 0.01
            timeend = 0.5

            @info (ArrayType, FT, dim, linearsolvertype, l)
            vtkdir = output ? "vtk_advection" *
                              "_poly$(polynomialorder)" *
                              "_dim$(dim)_$(ArrayType)_$(FT)" *
                              "_$(linearsolvertype)_level$(l)" : nothing
            result[l] = run(mpicomm, ArrayType, dim, topl, polynomialorder,
                            timeend, FT, dt, n, α, β, μ, δ, vtkdir,
                            outputtime, linearsolvertype)
            @test result[l] ≈ FT(expected_result[dim, l, FT])
          end
          @info begin
            msg = ""
            for l = 1:numlevels-1
              rate = log2(result[l]) - log2(result[l+1])
              msg *= @sprintf("\n  rate for level %d = %e\n", l, rate)
            end
            msg
          end
        end
      end
    end
  end
end

nothing

