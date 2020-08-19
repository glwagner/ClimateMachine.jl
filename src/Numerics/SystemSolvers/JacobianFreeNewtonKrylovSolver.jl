using Printf

export BatchedJacobianFreeNewtonKrylovSolver

mutable struct JacobianAction{F, FT}
    f!::F
    ϵ::FT
    cache_Fqdq
    cache_Fq
end

"""
Approximations the action of the Jacobian of a nonlinear
form on a vector `Δq` using the difference quotient:

∂F(q)      F(q + ϵΔq) - F(q)
---- Δq ≈ -------------------
 ∂q                ϵ

"""

function (op::JacobianAction)(dQ, Q, args...)
    f! = op.f!
    Fq = op.cache_Fq
    Fqdq = op.cache_Fqdq
    ϵ = op.ϵ

    f!(Fq, Q, args..., increment = false)
    f!(Fqdq, Q + ϵ .* dQ, args..., increment = false)
    dQ .= (Fqdq .- Fq) ./ ϵ
end

mutable struct BatchedJacobianFreeNewtonKrylovSolver{ET, TT, AT} <: AbstractNonlinearSolver
    # Tolerances
    ϵ::ET
    tol::TT
    # Max number of Newton iterations
    M::Int
    # Linear solver for the Jacobian system
    linearsolver
    # residual
    residual::AT
end

function BatchedJacobianFreeNewtonKrylovSolver(
    Q,
    linearsolver;
    ϵ = 1.e-8,
    tol = 1.e-6,
    M = 30,
)
    FT = eltype(Q)
    @info "Before"
    residual = similar(Q)
    @info "After"
    return BatchedJacobianFreeNewtonKrylovSolver(FT(ϵ), FT(tol), M, linearsolver, residual)
end

function initialize!(
    rhs!,
    Q,
    Qrhs,
    solver::BatchedJacobianFreeNewtonKrylovSolver,
    args...,
)
    # where R = Qrhs - F(Q)
    R = solver.residual
    # Computes F(Q) and stores in R
    rhs!(R, Q, args...)
    # Computes R = R - Qrhs
    R .-= Qrhs
    return norm(R, weighted_norm)
end

# rhs!(Q) =  F(Q)
# jvp!(Q)  = J(Q)ΔQ, here J(Q) = dF
# factors(Q) is the approximation of J(Q)
function donewtoniteration!(
    rhs!,   
    jvp!,                
    factors,
    Q,
    Qrhs,
    solver::BatchedJacobianFreeNewtonKrylovSolver,
    args...,
)

    FT = eltype(Q)
    ΔQ = similar(Q)
    ΔQ .= FT(0.0)

    #############################################################
    # Groupsize = "number of threads"
    # WANT: Execute N independent Newton iterations, where
    # N = Groupsize.
    #=

    R(Q) == 0, R = N - Qrhs, where N = rhs!

    N(Q) = Q - V(Q), where V(Q) is the 1-D nonlinear operator

    =#

    # Compute right-hand side for Jacobian system:
    # J(Q)ΔQ = -R
    # where R = Qrhs - F(Q)
    R = solver.residual
    # Computes F(Q) and stores in R
    rhs!(R, Q, args...)
    # Computes R = R - Qrhs
    R .-= Qrhs
    r0norm = norm(R, weighted_norm)
    # @info "Initial nonlinear residual F(Q): $r0norm"

    #= Inside linearsolve!
        need to perform the following:
        1. Preconditioning step (right PC)
            - Solve for w: Pw = ΔQ, where P
            is our preconditioner
        2. apply jvp on w
    =#

    # factors is an approximation of J(Q)

    iters = linearsolve!(
        jvp!,
        factors,
        solver.linearsolver,
        ΔQ,
        -R,
        args...,
    )

    # Newton correction
    Q .+= ΔQ

    # Reevaluate residual with new solution
    rhs!(R, Q, args...)
    R .-= Qrhs
    resnorm = norm(R, weighted_norm)
    # @info "Nonlinear residual F(Q) after solving Jacobian system: $resnorm"
    #############################################################
    
    return resnorm, iters
end
