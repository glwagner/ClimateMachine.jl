# Given an expression that is either a symbol which is a type name or a
# Union of types, each of which shares the same supertype, return the
# supertype.
parent_type(sym::Symbol) = supertype(getfield(@__MODULE__, sym))
function parent_type(ex::Expr)
    @assert ex.head == :curly && ex.args[1] == :Union
    st = supertype(getfield(@__MODULE__, ex.args[2]))
    for ut in ex.args[3:end]
        otherst = supertype(getfield(@__MODULE__, ut))
        @assert otherst == st
    end
    return st
end
# Return `true` if the specified symbol is a type name that is a subtype
# of `BalanceLaw` and `false` otherwise.
isa_bl(sym::Symbol) = any(
    bl -> endswith(bl, "." * String(sym)),
    map(bl -> String(Symbol(bl)), subtypes(BalanceLaw)),
)
isa_bl(ex) = false

# Generate a `VariableTemplates` defining function.
function generate_vars_fun(varsfun, dn, DT, vardecls, compdecls)
    quote
        function $varsfun($dn::$DT, FT)
            @vars begin
                $(vardecls...)
                $(compdecls...)
            end
        end
    end
end

# Generate the `VariableTemplates` defining functions for the specified
# `dvarnames`.
function generate_vars_funs(name, config_type, on_grid, params_type, dvarnames)
    varsfun = Symbol("vars_", name)
    CT = getfield(ConfigTypes, config_type)

    # Group the diagnostic variables specified by the dispatch types
    # specified in their implementations.
    DT_name_map = Dict()
    DT_var_map = Dict()
    for dvname in dvarnames
        dvar = AllDiagnosticVars[CT][dvname]
        dispatch_arg = first(dv_args(CT(), dvar))
        DT = dispatch_arg[2]
        compname = get(DT_name_map, DT, nothing)
        if isnothing(compname)
            DT_name_map[DT] = dispatch_arg[1]
        else
            @assert compname == dispatch_arg[1]
        end
        dvlst = get!(DT_var_map, DT, [])
        push!(dvlst, :($(Symbol(dvname))::FT))
    end

    # Generate a function for each dispatch type.
    vars_funs = []
    for (DT, dvlst) in DT_var_map
        # The top-level function must be for the `BalanceLaw` and it must
        # also contain component declarations, i.e. calls to the other
        # generated functions.
        complst = []
        if isa_bl(DT)
            for (otherDT, compname) in DT_name_map
                if otherDT != DT
                    push!(
                        complst,
                        :($(compname)::$varsfun(bl.$compname, FT)),
                    )
                end
            end
        else
            # Add an empty function for the parent of the dispatch type.
            push!(
                vars_funs,
                generate_vars_fun(
                    varsfun,
                    DT_name_map[DT],
                    parent_type(DT),
                    [],
                    [],
                ),
            )
        end
        push!(
            vars_funs,
            generate_vars_fun(
                varsfun,
                DT_name_map[DT],
                DT,
                dvlst,
                complst,
            ),
        )
    end

    return Expr(:block, (vars_funs...))
end

# Generate the `dims` dictionary for `Writers.init_data`.
function generate_init_dims(name, config_type, on_grid, dvarnames, dvars)
    # Set up an error for when no InterpolationTopology is specified but the
    # group is on an interpolated grid.
    err_ex = quote end
    if on_grid === :GridInterpolated
        err_ex = quote
            throw(
                ArgumentError(
                    "$name is on GridInterpolated and requires " *
                    "an InterpolationTopology"
                )
            )
        end
    end

    # For a diagnostics group run on the DG grid, we add some
    # "dimensions". For pointwise diagnostics, we add `nodes` and
    # `elements`. For horizontal averages, we add a `z` dimension.
    add_dim_ex = quote end
    if on_grid === :GridDG
        add_z_dim_ex = quote end
        if any(dv -> typeof(dv) <: HorizontalAverage, dvars)
            add_z_dim_ex = quote
                #TODO uncomment dims["z"] = (AtmosCollected.zvals, Dict())
                dims["z"] = (collect(1:100), Dict())
            end
        end
        add_ne_dims_ex = quote end
        if any(dv -> typeof(dv) <: PointwiseDiagnostic, dvars)
            add_ne_dims_ex = quote
                dims["nodes"] = (collect(1:npoints), Dict())
                dims["elements"] = (collect(1:nrealelem), Dict())
            end
        end
        add_dim_ex = quote
            $(add_z_dim_ex)
            $(add_ne_dims_ex)
        end
    end

    quote
        dims = dimensions(interpol)
        if isempty(dims)
            $(err_ex)
            $(add_dim_ex)
        elseif interpol isa InterpolationCubedSphere
            # Adjust `level` on the sphere.
            level_val = dims["level"]
            dims["level"] = (
                level_val[1] .- FT(planet_radius(Settings.param_set)),
                level_val[2],
            )
        end
        dims
    end
end

# Generate the `vars` dictionary for `Writers.init_data`.
function generate_init_vars(name, config_type, on_grid, dvarnames, dvars)
    CT = getfield(ConfigTypes, config_type)

    elems = Any[]
    for dvar in dvars
        rhs = :((dv_dims($dvar, dims), FT, $(dv_attrib(CT(), dvar))))
        lhs = :($(dv_name(CT(), dvar)))
        push!(elems, :($lhs => $rhs))
    end

    quote
        OrderedDict($(Expr(:tuple, elems...))...)
    end
end

# Generate `Diagnostics.$(name)_init(...)` which will initialize the
# `DiagnosticsGroup` when called.
function generate_init_fun(name, config_type, on_grid, params_type, dvarnames)
    init_name = Symbol(name, "_init")
    CT = getfield(ConfigTypes, config_type)

    dvars = [AllDiagnosticVars[CT][dvname] for dvname in dvarnames]

    quote
        function $init_name(dgngrp, curr_time)
            interpol = dgngrp.interpol
            mpicomm = Settings.mpicomm
            mpirank = 0 #MPI.Comm_rank(mpicomm)
            dg = Settings.dg
            bl = dg.balance_law
            grid = dg.grid
            topology = grid.topology
            N = polynomialorder(grid)
            Nq = N + 1
            Nqk = dimensionality(grid) == 2 ? 1 : Nq
            npoints = Nq * Nq * Nqk
            nrealelem = length(topology.realelems)
            nvertelem = topology.stacksize
            nhorzelem = div(nrealelem, nvertelem)
            Q = Settings.Q
            FT = eltype(Q)

            if dgngrp.onetime
                atmos_collect_onetime(Settings.mpicomm, Settings.dg, Settings.Q)
            end

            if mpirank == 0
                dims = $(generate_init_dims(name, config_type, on_grid, dvarnames, dvars))
                vars = $(generate_init_vars(name, config_type, on_grid, dvarnames, dvars))
                println(dims)
                println(vars)
                # create the output file
                dprefix = @sprintf(
                    "%s_%s_%s",
                    dgngrp.out_prefix,
                    dgngrp.name,
                    Settings.starttime,
                )
                dfilename = joinpath(Settings.output_dir, dprefix)
                #init_data(dgngrp.writer, dfilename, dims, vars)
            end

            return nothing
        end
    end
end

# Generate `Diagnostics.$(name)_collect(...)` which when called,
# performs a collection of all the diagnostic variables in the group
# and writes them out.
function generate_collect_fun_dev(name, config_type, on_grid, params_type, dvarnames)
    collect_name = Symbol(name, "_collect")
    quote
        function $collect_name(dgngrp, curr_time)
            mpicomm = Settings.mpicomm
            mpirank = MPI.Comm_rank(mpicomm)
            dg = Settings.dg
            bl = dg.balance_law
            Q = Settings.Q
            FT = eltype(Q)
            interpol = dgngrp.interpol

            intermediates =
                $(generate_intermediates(name, config_type, dvarnames))

            $(generate_collect_vars(name, ))
            if mpirank == 0
                varvals = OrderedDict()
                # XXX
                append_data(dgngrp.writer, varvals, curr_time)
            end

            return nothing
        end
    end
end
function generate_collect_fun(name, config_type, on_grid, params_type, dvarnames)
    collect_name = Symbol(name, "_collect")
    quote
        function $collect_name(dgngrp, curr_time)
            mpicomm = Settings.mpicomm
            mpirank = MPI.Comm_rank(mpicomm)
            dg = Settings.dg
            bl = dg.balance_law
            Q = Settings.Q
            FT = eltype(Q)
            interpol = dgngrp.interpol

            if mpirank == 0
                varvals = OrderedDict()
                # XXX
                append_data(dgngrp.writer, varvals, curr_time)
            end

            return nothing
        end
    end
end


# Generate `Diagnostics.$(name)_fini(...)`, which does nothing
# right now.
function generate_fini_fun(name, config_type, on_grid, params_type, dvarnames)
    fini_name = Symbol(name, "_fini")
    quote
        function $fini_name(dgngrp, curr_time) end
    end
end

# Generate `setup_$(name)(...)` which will create the `DiagnosticsGroup`
# for $name when called.
function generate_setup_fun(name, config_type, on_grid, params_type)
    setupfun = Symbol("setup_", name)
    initfun = Symbol(name, "_init")
    collectfun = Symbol(name, "_collect")
    finifun = Symbol(name, "_fini")
    quote
        function $setupfun(
            ::$config_type,
            params::$params_type,
            interval::String,
            out_prefix::String,
            writer = NetCDFWriter(),
            interpol = nothing,
        ) where {
            $config_type <: ClimateMachineConfigType,
            $params_type <: Union{Nothing, DiagnosticsGroupParams},
        }
            return DiagnosticsGroup(
                $(name),
                $(initfun),
                $(collectfun),
                $(finifun),
                interval,
                out_prefix,
                writer,
                interpol,
                $(on_grid === :GridDG),
                params,
            )
        end
    end
end
