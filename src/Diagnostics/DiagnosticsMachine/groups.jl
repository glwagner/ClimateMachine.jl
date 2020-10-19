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
mutable struct DiagnosticsGroup{DGP <: Union{Nothing, DiagnosticsGroupParams}}
    name::String
    init::Function
    collect::Function
    fini::Function
    interval::String
    out_prefix::String
    writer::AbstractWriter
    interpol::Union{Nothing, InterpolationTopology}
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
        params = nothing,
    ) = new{typeof(params)}(
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
    dvars...,
)
    setup = generate_setup(name, config_type, on_grid, params_type)
    dump(setup)
    init = generate_init(name, config_type, on_grid, params_type, dvars)
    dump(init)
    collect = generate_collect(name, config_type, on_grid, params_type, dvars)
    dump(collect)
    fini = generate_fini(name, dvars)
    dump(fini)

    return Expr(
        :block,
        esc(setup),
        esc(init),
        esc(collect),
        esc(fini),
    )
end

include("group_gen.jl")