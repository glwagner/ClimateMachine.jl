# Generate `setup_$(name)(...)` which will create the `DiagnosticsGroup`
# for $name when called.
function generate_setup(name, config_type, on_grid, params_type)
    setup_name = Symbol("setup_", name)
    init_name = Symbol(name, "_init")
    collect_name = Symbol(name, "_collect")
    fini_name = Symbol(name, "_fini")
    quote
        function $(setup_name)(
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
                Diagnostics.$(init_name),
                Diagnostics.$(collect_name),
                Diagnostics.$(fini_name),
                interval,
                out_prefix,
                writer,
                interpol,
                $(on_grid ∈ GridDG),
                params,
            )
        end
    end
end

# Generate the `dims` dictionary for `Writers.init_data`.
function generate_init_dims(name, config_type, on_grid, dvars)
    # Set up an error for when no InterpolationTopology is specified but the
    # group is on an interpolated grid.
    err_ex = quote end
    if on_grid ∈ GridInterpolated
        err_ex = quote
            throw(
                ArgumentError(
                    "$name is on GridInterpolated and requires " *
                    "an InterpolationTopology"
                )
            )
        end
    end

    # For a diagnostics group run on the DG grid, we add some "dimensions".
    # For horizontal averages, we add a `z` dimension. For pointwise
    # diagnostics, we add `nodes` and `elements`.
    add_dim_ex = quote end
    if on_grid ∈ GridDG
        add_ne_dims_ex = quote end
        if any(dv -> dv isa PointwiseDiagnostic, dvars)
            add_ne_dims_ex = quote
                $(esc(dims))["nodes"] = (collect(1:$(esc(npoints))), Dict())
                $(esc(dims))["elements"] = (collect(1:$(esc(nrealelem))), Dict())
            end
        end
        add_z_dim_ex = quote end
        if any(dv -> dv isa HorizontalAverage, dvars)
            add_z_dim_ex = quote
                $(esc(dims))["z"] = ($(esc(AtmosCollected.zvals)), Dict())
            end
        end
        add_dim_ex = quote
            $(esc(dims)) = OrderedDict()
            $(add_ne_dims_ex)
            $(add_z_dim_ex)
        end
    end

    quote
        dims = dimensions($(esc(interpol)))
        if isempty(dims)
            $(err_ex)
            $(add_dim_ex)
        elseif $(esc(interpol)) isa InterpolationCubedSphere
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
function generate_init_vars(name, config_type, on_grid, dvars, dims)
    elems = ()

    for dvar in dvars
        elems = (
            elems...,
            dv_name(config_type, dvar) => (
                :(dv_dims($dvar, dims)),
                FT,
                dv_attrib(config_type, dvar),
            )
        )
    end

    quote
        OrderedDict($(elems...))
    end
end

# Generate `Diagnostics.$(name)_init(...)` which will initialize the
# `DiagnosticsGroup` when called.
function generate_init(name, config_type, on_grid, params_type, dvars)
    init_name = Symbol(name, "_init")
    quote
        function $(esc(init_name))(dgngrp, curr_time)
            mpicomm = Settings.mpicomm
            mpirank = MPI.Comm_rank(mpicomm)
            dg = Settings.dg
            bl = dg.balance_law
            Q = Settings.Q
            FT = eltype(Q)
            interpol = dgngrp.interpol

            if dgngrp.onetime
                atmos_collect_onetime(Settings.mpicomm, Settings.dg, Settings.Q)
            end

            if mpirank == 0
                dims = $(generate_init_dims(name, config_type, on_grid, dvars))
                vars = $(generate_init_vars(name, config_type, on_grid, dvars))

                # create the output file
                dprefix = @sprintf(
                    "%s_%s_%s",
                    dgngrp.out_prefix,
                    dgngrp.name,
                    Settings.starttime,
                )
                dfilename = joinpath(Settings.output_dir, dprefix)
                init_data(dgngrp.writer, dfilename, dims, vars)
            end

            return nothing
        end
    end
end

# Generate `Diagnostics.$(name)_collect(...)` which when called,
# performs a collection of all the diagnostic variables in the group
# and writes them out.
function generate_collect(name, dvars)
    collect_name = Symbol(name, "_collect")
    quote
        function $(esc(collect_name))(dgngrp, curr_time)
            mpicomm = Settings.mpicomm
            mpirank = MPI.Comm_rank(mpicomm)
            dg = Settings.dg
            bl = dg.balance_law
            Q = Settings.Q
            FT = eltype(Q)
            interpol = dgngrp.interpol

            intermediates =
                $(generate_intermediates(name, config_type, dvars))

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

# Generate `Diagnostics.$(name)_fini(...)`, which does nothing
# right now.
function generate_fini(name, vars)
    fini_name = Symbol(name, "_fini")
    quote
        function $(esc(fini_name))(dgngrp, curr_time) end
    end
end
