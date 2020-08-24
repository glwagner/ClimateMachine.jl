using Printf

export BatchedJacobianFreeNewtonKrylovSolver, JacobianAction

mutable struct JacobianAction{FT, AT}
    rhs!
    ϵ::FT
    Q::AT
    cache_Fq::AT
    cache_Fqdq::AT
end

function JacobianAction(rhs!, Q, ϵ)
    cache_Fq = similar(Q)
    cache_Fqdq = similar(Q)
    return JacobianAction(rhs!, ϵ, Q, cache_Fq, cache_Fqdq)
end


"""
Approximations the action of the Jacobian of a nonlinear
form on a vector `Δq` using the difference quotient:

∂F(q)      F(q + ϵΔq) - F(q)
---- Δq ≈ -------------------
 ∂q                ϵ

"""

function (op::JacobianAction)(JΔQ, dQ, args...)
    rhs! = op.rhs!
    Q = op.Q
    ϵ = op.ϵ
    Fq = op.cache_Fq
    Fqdq = op.cache_Fqdq
    
    FT = eltype(dQ)
    n = length(dQ)
    normdQ = norm(dQ, weighted_norm)

    if normdQ > ϵ
        factor = FT(1 / (n*normdQ))
    else
        # initial newton step, ΔQ = 0
        factor = FT(1 / n)
    end

    β = √ϵ
    e = factor * β * sum(abs.(Q)) + β

    rhs!(Fqdq, Q .+ e .* dQ, args...)

    JΔQ .= (Fqdq .- Fq) ./ e
end

function update_Q!(op::JacobianAction, Q, args...)
    op.Q .= Q
    Fq = op.cache_Fq

    op.rhs!(Fq, Q, args...)
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
    residual = similar(Q)
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
