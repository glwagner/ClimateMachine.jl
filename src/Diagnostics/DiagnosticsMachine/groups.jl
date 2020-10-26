"""
    DiagnosticsGroupParams

Base type for any extra paramters needed by a diagnostics group.
"""
abstract type DiagnosticsGroupParams end

"""
    DiagnosticsGroup

A diagnostics group contains a set of diagnostic variables that are
collected at the same interval and written to the same file. A group's
variables are collected either on the DG grid or on an interpolated
grid.
"""
mutable struct DiagnosticsGroup{
    DGP <: Union{Nothing, DiagnosticsGroupParams},
    DGI <: Union{Nothing, InterpolationTopology},
}
    name::String
    init::Function
    collect::Function
    fini::Function
    interval::String
    out_prefix::String
    writer::AbstractWriter
    interpol::DGI
    onetime::Bool
    params::DGP

    DiagnosticsGroup(
        name,
        init,
        collect,
        fini,
        interval,
        out_prefix,
        writer,
        interpol,
        onetime,
        params,
    ) = new{typeof(params), typeof(interpol)}(
        name,
        init,
        collect,
        fini,
        interval,
        out_prefix,
        writer,
        interpol,
        onetime,
        params,
    )
end

# `GenericCallbacks` implementations for `DiagnosticsGroup`s
function GenericCallbacks.init!(dgngrp::DiagnosticsGroup, solver, Q, param, t)
    @info @sprintf(
        """
    Diagnostics: %s
        initializing at %8.2f""",
        dgngrp.name,
        t,
    )
    dgngrp.init(dgngrp, t)
    dgngrp.collect(dgngrp, t)
    return nothing
end
function GenericCallbacks.call!(dgngrp::DiagnosticsGroup, solver, Q, param, t)
    @tic diagnostics
    @info @sprintf(
        """
    Diagnostics: %s
        collecting at %8.2f""",
        dgngrp.name,
        t,
    )
    dgngrp.collect(dgngrp, t)
    @toc diagnostics
    return nothing
end
function GenericCallbacks.fini!(dgngrp::DiagnosticsGroup, solver, Q, param, t)
    @info @sprintf(
        """
    Diagnostics: %s
        finishing at %8.2f""",
        dgngrp.name,
        t,
    )
    dgngrp.collect(dgngrp, t)
    dgngrp.fini(dgngrp, t)
    return nothing
end

"""
    OnGridKind

The variables in a diagnostic group are computed from either the DG
grid or an interpolated grid.
"""
abstract type OnGridKind end
struct GridDG <: OnGridKind end
struct GridInterpolated <: OnGridKind end

"""
    @diagnostics_group

Generate the functions needed to establish and use a `DiagnosticsGroup`
containing the named `DiagnosticVar`s.
"""
macro diagnostics_group(
    name,
    config_type,
    on_grid,
    params_type,
    dvarnames...,
)
    ex = generate_vars_funs(name, config_type, on_grid, params_type, dvarnames)
    vars_funs = esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, ex)))
    #println(vars_funs)
    init_fun = esc(generate_init_fun(name, config_type, on_grid, params_type, dvarnames))
    #init_fun = esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, ex)))
    println(init_fun)
    ex = generate_collect_fun(name, config_type, on_grid, params_type, dvarnames)
    collect_fun = esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, ex)))
    #println(collect_fun)
    ex = generate_fini_fun(name, config_type, on_grid, params_type, dvarnames)
    fini_fun = esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, ex)))
    #println(fini_fun)
    ex = generate_setup_fun(name, config_type, on_grid, params_type)
    setup_fun = esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, ex)))
    #println(setup_fun)

    return Expr(
        :block,
        vars_funs,
        setup_fun,
        init_fun,
        collect_fun,
        fini_fun,
    )
end

include("group_gen.jl")