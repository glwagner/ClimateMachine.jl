export LinearBackwardEulerSolver, AbstractBackwardEulerSolver
export NonLinearBackwardEulerSolver

abstract type AbstractImplicitOperator end

"""
    op! = EulerOperator(f!, ϵ)

Construct a linear operator which performs an explicit Euler step ``Q + α
f(Q)``, where `f!` and `op!` both operate inplace, with extra arguments passed
through, i.e.
```
op!(LQ, Q, args...)
```
is equivalent to
```
f!(dQ, Q, args...)
LQ .= Q .+ ϵ .* dQ
```
"""
mutable struct EulerOperator{F, FT} <: AbstractImplicitOperator
    f!::F
    ϵ::FT
end

function (op::EulerOperator)(LQ, Q, args...)
    op.f!(LQ, Q, args..., increment = false)
    @. LQ = Q + op.ϵ * LQ
end

"""
    AbstractBackwardEulerSolver

An abstract backward Euler method
"""
abstract type AbstractBackwardEulerSolver end

"""
    (be::AbstractBackwardEulerSolver)(Q, Qhat, α, param, time)

Each concrete implementations of `AbstractBackwardEulerSolver` should provide a
callable version which solves the following system for `Q`
```
    Q = Qhat + α f(Q, param, time)
```
where `f` is the ODE tendency function, `param` are the ODE parameters, and
`time` is the current ODE time. The arguments `Q` should be modified in place
and should not be assumed to be initialized to any value.
"""
(be::AbstractBackwardEulerSolver)(Q, Qhat, α, p, t) =
    throw(MethodError(be, (Q, Qhat, α, p, t)))

"""
    Δt_is_adjustable(::AbstractBackwardEulerSolver)

Return `Bool` for whether this backward Euler solver can be updated. default is
`false`.
"""
Δt_is_adjustable(::AbstractBackwardEulerSolver) = false

"""
    update_backward_Euler_solver!(::AbstractBackwardEulerSolver, α)

Update the given backward Euler solver for the parameter `α`; see
['AbstractBackwardEulerSolver'](@ref). Default behavior is no change to the
solver.
"""
update_backward_Euler_solver!(::AbstractBackwardEulerSolver, Q, α) = nothing

"""
    setup_backward_Euler_solver(solver, Q, α, tendency!)

Returns a concrete implementation of an `AbstractBackwardEulerSolver` that will
solve for `Q` in systems of the form of
```
    Q = Qhat + α f(Q, param, time)
```
where `tendency!` is the in-place tendency function. Not the array `Q` is just
passed in for type information, e.g., `Q` the same `Q` will not be used for all
calls to the solver.
"""
setup_backward_Euler_solver(solver::AbstractBackwardEulerSolver, _...) = solver

"""
    LinearBackwardEulerSolver(::AbstractSystemSolver; isadjustable = false)

Helper type for specifying building a backward Euler solver with a linear
solver.  If `isadjustable == true` then the solver can be updated with a new
time step size.
"""
struct LinearBackwardEulerSolver{LS}
    solver::LS
    isadjustable::Bool
    preconditioner::Bool
    LinearBackwardEulerSolver(solver; isadjustable = false, preconditioner = false) =
        new{typeof(solver)}(solver, isadjustable, preconditioner)
end

"""
    LinBESolver

Concrete implementation of an `AbstractBackwardEulerSolver` to use linear
solvers of type `AbstractSystemSolver`. See helper type
[`LinearBackwardEulerSolver`](@ref)
```
    Q = Qhat + α f(Q, param, time)
```
"""
mutable struct LinBESolver{FT, F, LS} <: AbstractBackwardEulerSolver
    α::FT
    f_imp!::F
    solver::LS
    isadjustable::Bool
    preconditioner::Bool
    factors
end

Δt_is_adjustable(lin::LinBESolver) = lin.isadjustable

function setup_backward_Euler_solver(lin::LinearBackwardEulerSolver, Q, α, f_imp!)
    FT = eltype(α)
    rhs! = EulerOperator(f_imp!, -α)

    factors = prefactorize(rhs!, lin.solver, Q, nothing, FT(NaN))

    LinBESolver(α, f_imp!, lin.solver, lin.isadjustable, lin.preconditioner, factors)
end

function update_backward_Euler_solver!(lin::LinBESolver, Q, α)
    lin.α = α
    FT = eltype(Q)
    # for direct solver, update factors
    # for iterative solver, set factors to Nothing (TODO optimize)
    lin.factors = prefactorize(EulerOperator(lin.f_imp!, -α), lin.solver, Q, nothing, FT(NaN),)
end

function (lin::LinBESolver)(Q, Qhat, α, p, t)
    rhs! = EulerOperator(lin.f_imp!, -α)

    if lin.α != α
        @assert lin.isadjustable
        update_backward_Euler_solver!(lin, Q, α)
    end

    if lin.preconditioner
        # update the preconditioner, lin.factors
        FT = eltype(α)
        # TODO what is the single_column
        single_column = false
        lin.factors = preconditioner(rhs!, single_column, Q, nothing, FT(NaN), )
    end
    linearsolve!(rhs!, lin.factors, lin.solver, Q, Qhat, p, t)
end


###################################################################################################
struct NonLinearBackwardEulerSolver{NLS}
    nlsolver::NLS
    isadjustable::Bool
    # preconditioner
    preconditioner::Bool
    # linearized operator, does not have α information
    lin_f_imp!
    function NonLinearBackwardEulerSolver(nlsolver; isadjustable = false, preconditioner = false, lin_f_imp! = nothing)
        NLS = typeof(nlsolver)
        return new{NLS}(nlsolver, isadjustable, preconditioner, lin_f_imp!)
    end
end

mutable struct NonLinBESolver{FT, F, NLS} <: AbstractBackwardEulerSolver
    α::FT
    f_imp!::F
    nlsolver::NLS
    isadjustable::Bool
    # preconditioner
    preconditioner::Bool
    # linearized operator, does not have α information
    lin_f_imp!
    # preconditioner factors, has α information
    factors
    
end

Δt_is_adjustable(nlsolver::NonLinBESolver) = nlsolver.isadjustable

function setup_backward_Euler_solver(
    nlsolver::NonLinearBackwardEulerSolver,
    Q,
    α,
    f_imp!,
)   
    @info "NonLinearBackwardEulerSolver setup_backward_Euler_solver : ",  α
    FT = eltype(α)
    NonLinBESolver(
        α,
        f_imp!,
        nlsolver.nlsolver,
        nlsolver.isadjustable,
        nlsolver.preconditioner, 
        nlsolver.lin_f_imp!,
        nothing
    )
end

# function update_backward_Euler_solver!(nlbesolver::NonLinBESolver, Q, α)
#     # TODO: What else needs to be updated? The linear solver?
#     # Should the `NonLinBESolver` object also point to
#     # a corresponding `LinBESolver` for the linear problem?

#     @info "NonLinBESolver update_backward_Euler_solver!"
#     nlbesolver.α = α
#     update_backward_Euler_solver!(nlbesolver.nlsolver.lsolver, Q, α)
# end

function (nlbesolver::NonLinBESolver)(Q, Qhat, α, p, t)

    # if nlbesolver.α != α

    #     @info "nlbesolver.α != α"
    #     @assert nlbesolver.isadjustable
    #     update_backward_Euler_solver!(nlbesolver, Q, α)
    # end

    # rhs! = EulerOperator(nlbesolver.f_imp!, -α)

    rhs! = nlbesolver.f_imp!

    # if  nlbesolver.preconditioner
    #     # update the preconditioner, lin.factors
    #     FT = eltype(α)
    #     # TODO what is the single_column
    #     single_column = false
    #     nlbesolver.factors = preconditioner(EulerOperator(nlbesolver.lin_f_imp!, -α), single_column, Q, nothing, FT(NaN), )
    # end

    # Call "solve" function in SystemSolvers
    nonlinearsolve!(
        rhs!,
        nlbesolver.nlsolver,
        Q,
        Qhat,
        p,
        t,
    )
end
