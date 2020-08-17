module SystemSolvers

using ..MPIStateArrays
using ..MPIStateArrays: array_device, realview

using ..Mesh.Grids
import ..Mesh.Grids: polynomialorder, dimensionality
using ..Mesh.Topologies
using ..DGMethods
using ..DGMethods: DGModel
using ..BalanceLaws

using Adapt
using LinearAlgebra
using LazyArrays
using StaticArrays
using KernelAbstractions

const weighted_norm = false

# just for testing SystemSolvers
LinearAlgebra.norm(A::MVector, p::Real, weighted::Bool) = norm(A, p)
LinearAlgebra.norm(A::MVector, weighted::Bool) = norm(A, 2, weighted)
LinearAlgebra.dot(A::MVector, B::MVector, weighted) = dot(A, B)
LinearAlgebra.norm(A::AbstractVector, p::Real, weighted::Bool) = norm(A, p)
LinearAlgebra.norm(A::AbstractVector, weighted::Bool) = norm(A, 2, weighted)
LinearAlgebra.dot(A::AbstractVector, B::AbstractVector, weighted) = dot(A, B)

export linearsolve!, settolerance!, prefactorize, preconditioner, preconditioner_solve!, preconditioner_matprod!
export AbstractSystemSolver, AbstractIterativeSystemSolver
export nonlinearsolve!

"""
    AbstractSystemSolver

This is an abstract type representing a generic linear solver.
"""
abstract type AbstractSystemSolver end

abstract type AbstractNonlinearSolver <: AbstractSystemSolver end

struct LSOnly <: AbstractNonlinearSolver
    linearsolver
end

function donewtoniteration!(linearoperator!, Q, Qrhs, solver::LSOnly, tol, args...)
    linearsolve!(
        linearoperator!,
        nothing,
        solver.linearsolver,
        Q,
        Qrhs,
        args...;
        max_iters = getmaxiterations(solver.linearsolver),
    )
end

function apply_jacobian!(
    JΔQ,
    implicitoperator!,
    Q,
    dQ,
    ϵ,
    args...,
)
    n = length(dQ)
    normdQ = norm(dQ, weighted_norm)

    if normdQ > ϵ
        factor = (1 / (n*normdQ))
    else
        # initial newton step, ΔQ = 0
        factor = 1 / n
    end

    β = √ϵ
    e = factor * β * sum(abs.(Q)) + β

    Fq = similar(Q)
    Fqdq = similar(Q)
    implicitoperator!(Fq, Q, args..., increment = false)
    implicitoperator!(Fqdq, Q .+ e .* dQ, args..., increment = false)
    JΔQ .= (Fqdq .- Fq) ./ e
end

"""

Solving F(Q) == 0 via Newton,

where `F = N(Q) - Qrhs`, N(Q) is
`implicitoperator!`.

"""
function nonlinearsolve!(
    solver::AbstractNonlinearSolver,
    Q::AT,
    Qrhs,
    args...;
    max_newton_iters = 10,
    cvg = Ref{Bool}(),
) where {AT}

    tol = solver.tol
    converged = false
    iters = 0
    implicitoperator! = solver.rhs!

    # Initialize NLSolver, compute initial residual
    initial_residual_norm =
        initialize!(implicitoperator!, Q, Qrhs, solver, args...)
    if initial_residual_norm < tol
        converged = true
    end
    converged && return iters

    """
    Want:
        (*) linearoperator!(Result, CurrentState:ΔQ, args...)
    
    Want the Jacobian action (jvp!) to behave just like
    a standard rhs evaluation as in (*)
    """
    # Create Jacobian action here?
    jvp! = (JΔQ, ΔQ, args...) -> apply_jacobian!(JΔQ, 
        implicitoperator!,
        Q,
        ΔQ,
        solver.ϵ,
        args...,
    )

    while !converged && iters < max_newton_iters
        residual_norm, linear_iterations =
            donewtoniteration!(implicitoperator!, jvp!, Q, Qrhs, solver, tol, args...)
        @info "Linear solver converged in $linear_iterations iterations"
        iters += 1

        if !isfinite(residual_norm)
            error("norm of residual is not finite after $iters iterations of `donewtoniteration!`")
        end

        # Check residual_norm / norm(R0)
        # Comment: Should we check "correction" magitude?
        # ||Delta Q|| / ||Q|| ?
        relresidual = residual_norm / initial_residual_norm
        # @info "Relative residual, current residual, initial residual: $relresidual, $residual_norm, $initial_residual_norm" 
        if relresidual < tol || residual_norm < tol
            @info "Newton converged in $iters iterations!"
            converged = true
        end
    end

    converged || @warn "Nonlinear solver did not attain convergence after $iters iterations"
    cvg[] = converged

    iters
end

"""
    AbstractIterativeSystemSolver

This is an abstract type representing a generic iterative
linear solver.

The available concrete implementations are:

  - [`GeneralizedConjugateResidual`](@ref)
  - [`GeneralizedMinimalResidual`](@ref)
"""
abstract type AbstractIterativeSystemSolver <: AbstractSystemSolver end

"""
    settolerance!(solver::AbstractIterativeSystemSolver, tolerance, relative)

Sets the relative or absolute tolerance of the iterative linear solver
`solver` to `tolerance`.
"""
settolerance!(
    solver::AbstractIterativeSystemSolver,
    tolerance,
    relative = true,
) = (relative ? (solver.rtol = tolerance) : (solver.atol = tolerance))

doiteration!(
    linearoperator!,
    factors,
    Q,
    Qrhs,
    solver::AbstractIterativeSystemSolver,
    tolerance,
    args...,
) = throw(MethodError(
    doiteration!,
    (linearoperator!, factors, Q, Qrhs, solver, tolerance, args...),
))

initialize!(
    linearoperator!,
    factors, 
    Q,
    Qrhs,
    solver::AbstractIterativeSystemSolver,
    args...,
) = throw(MethodError(initialize!, (linearoperator!, factors, Q, Qrhs, solver, args...)))

"""
    prefactorize(linop!, linearsolver, args...)

Prefactorize the in-place linear operator `linop!` for use with `linearsolver`.
"""
prefactorize(linop!, linearsolver::AbstractIterativeSystemSolver, args...) =
    linop!

"""
    linearsolve!(linearoperator!, solver::AbstractIterativeSystemSolver, Q, Qrhs, args...)

Solves a linear problem defined by the `linearoperator!` function and the state
`Qrhs`, i.e,

```math
L(Q) = Q_{rhs}
```

using the `solver` and the initial guess `Q`. After the call `Q` contains the
solution.  The arguments `args` is passed to `linearoperator!` when it is
called.

jvp! = (ΔQ) -> F(Q + ΔQ)
"""


function linearsolve!(
    linearoperator!,
    factors,
    solver::AbstractIterativeSystemSolver,
    Q,
    Qrhs,
    args...;
    max_iters = length(Q),
    cvg = Ref{Bool}(),
)
    converged = false
    iters = 0
    @info "I am in linearsolve!"
    converged, residual_norm0 = initialize!(linearoperator!, factors, Q, Qrhs, solver, args...)
    converged && return iters

    while !converged && iters < max_iters
        converged, inner_iters, residual_norm =
            doiteration!(linearoperator!, factors, Q, Qrhs, solver, args...)

        iters += inner_iters

        if !isfinite(residual_norm)
            error("norm of residual is not finite after $iters iterations of `doiteration!`")
        end

        achieved_tolerance = residual_norm / residual_norm0 * solver.rtol
    end

    converged || @warn "Solver did not attain convergence after $iters iterations"
    cvg[] = converged

    iters
end

@kernel function linearcombination!(Q, cs, Xs, increment::Bool)
    i = @index(Global, Linear)
    if !increment
        @inbounds Q[i] = -zero(eltype(Q))
    end
    @inbounds for j in 1:length(cs)
        Q[i] += cs[j] * Xs[j][i]
    end
end

include("JacobianFreeNewtonKrylovSolver.jl")
include("generalized_minimal_residual_solver.jl")
include("generalized_conjugate_residual_solver.jl")
include("conjugate_gradient_solver.jl")
include("columnwise_lu_solver.jl")
include("batched_generalized_minimal_residual_solver.jl")

end
