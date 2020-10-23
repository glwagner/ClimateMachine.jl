# TODO:
# - scalars and others from #1596, especially core (conditional)
# - diagnostics from a specialized kernel, e.g. vorticity

"""
    States

Composite of the various states, used as a parameter to diagnostic
collection functions.
"""
struct States{PS,GFS,AS,TS}
    prognostic::PS
    gradient_flux::GFS
    auxiliary::AS
    thermodynamic::TS
end

"""
    AbstractIntermediates

Base type for a composite that will be generated.
"""
abstract type AbstractIntermediates end

"""
    DiagnosticVar

The base type for diagnostic variables.
"""
abstract type DiagnosticVar end

"""
    dv_name(::Type{CT}, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }
"""
function dv_name end

"""
    dv_attrib(::Type{CT}, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }
"""
function dv_attrib end

"""
    dv_args(::Type{CT}, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns a tuple of the arguments specified by the implementation of the
diagnostic variables.
"""
function dv_args end

"""
    dv_dims(::DVT, out_dims::ODT) where {
        DVT <: DiagnosticVar,
        ODT <: Union{Nothing, Tuple},
    }

The `out_dims` parameter may be `nothing`, or a tuple
with the names of the dimensions specified for the output
(by the `InterpolationTopology` for instance).
"""
function dv_dims end

# Default method for variable attributes.
dv_attrib(
    ::ClimateMachineConfigType,
    ::DiagnosticVar,
) = Dict()

# Default method for variable dimensions.
dv_dims(::DiagnosticVar, ::Nothing) = ("nodes", "elements")

# Generate a standardized type name from the diagnostic variable name.
function dv_type_name(dvtype, config_type, name)
    let uppers_in(s) =
        foldl((f, c) -> isuppercase(c) ? f * c : f, String(s), init = "")
        return uppers_in(config_type) * "_" * uppers_in(dvtype) * "_" * name
    end
end

# Generate the type and interface functions for a diagnostic variable.
function generate_dv_interface(
    dvtype,
    config_type,
    name,
    units = "",
    long_name = "",
    standard_name = "",
)
    dvtypname = Symbol(dv_type_name(dvtype, config_type, name))
    attrib_ex = quote end
    if any(a -> a != "", [units, long_name, standard_name])
        attrib_ex = quote
            dv_attrib(::Type{$config_type}, ::Type{$dvtypname}) =
                OrderedDict(
                    "units" => $units,
                    "long_name" => $long_name,
                    "standard_name" => $standard_name,
                )
        end
    end
    quote
        struct $dvtypname <: $dvtype end
        DiagnosticsMachine.AllDiagnosticVars[$config_type][$name] = $dvtypname
        dv_name(::Type{$config_type}, ::Type{$dvtypname}) = $name
        $(attrib_ex)
    end
end

# Helper to generate the implementation function for one or more
# diagnostic variables.
function generate_dv_function(
    dvtype,
    config_type,
    names,
    impl,
)
    dvfun = Symbol("dv_", dvtype)
    dvtypname_args = map(
        n -> :(Type{$n}),
        map(n -> Symbol(dv_type_name(dvtype, config_type, n)), names),
    )
    @capture(impl, ((args__,),) -> (body_)) ||
        @capture(impl, (args_) -> (body_)) ||
            error("Bad implementation for $(esc(names[1]))")
    split_fun_args = map(splitarg, args)
    fun_args = map(a -> :($(a[1])::$(a[2])), split_fun_args)
    quote
        function dv_args(
            ::Type{$config_type},
            ::Union{$(dvtypname_args...)},
        )
            $split_fun_args
        end
        function $dvfun(
            ::Type{$config_type},
            ::Union{$(dvtypname_args...)},
            $(fun_args...),
        )
            $body
        end
    end
end

# Interface to generate an implementation function for one or more
# diagnostic variables.
macro diagnostic_impl(
    impl,
    dvtype,
    config_type,
    names...,
)
    generate_dv_function(
        dvtype,
        config_type,
        names,
        impl,
    )
end

# Diagnostic variable types and interfaces to create diagnostic variables
# of these types.

"""
    IntermediateValue
"""
abstract type IntermediateValue <: DiagnosticVar end
dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_dims(::IntermediateValue, out_dims::Tuple) = out_dims

macro intermediate_value(config_type, name)
    iex = generate_dv_interface(:IntermediateValue, config_type, name)
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro intermediate_value(impl, config_type, name)
    iex = quote
        $(generate_dv_interface(:IntermediateValue, config_type, name))
        $(generate_dv_function(
            :IntermediateValue,
            config_type,
            [name],
            impl,
        ))
    end
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro intermediate_values(impl, config_type, names...)
    exprs = [
        generate_dv_interface(
            :IntermediateValue,
            config_type,
            name,
        )
        for name in names
    ]
    fex = generate_dv_function(
        :IntermediateValue,
        config_type,
        names,
        impl,
    )
    push!(exprs, fex)
    iex = quote
        $(exprs...)
    end
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

"""
    PointwiseDiagnostic

A diagnostic with the same dimensions as the `from_grid` chosen in
the diagnostics group.
"""
abstract type PointwiseDiagnostic <: DiagnosticVar end
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_dims(::PointwiseDiagnostic, out_dims::Tuple) = out_dims

macro pointwise_diagnostic(config_type, name)
    iex = generate_dv_interface(PointwiseDiagnostic, config_type, name)
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro pointwise_diagnostic(
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = generate_dv_interface(
        :PointwiseDiagnostic,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro pointwise_diagnostic(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = quote
        $(generate_dv_interface(
            :PointwiseDiagnostic,
            config_type,
            name,
            units,
            long_name,
            standard_name,
        ))
        $(generate_dv_function(
            :PointwiseDiagnostic,
            config_type,
            [name],
            impl,
        ))
    end
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro pointwise_diagnostic_impl(
    impl,
    config_type,
    names...,
)
    iex = quote
        $(generate_dv_function(
            :PointwiseDiagnostic,
            config_type,
            names,
            impl,
        ))
    end
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

"""
    HorizontalAverage

A horizontal reduction into a single vertical dimension.
"""
abstract type HorizontalAverage <: DiagnosticVar end
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{HorizontalAverage}},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_dims(::HorizontalAverage, ::Nothing) = ("z",)
dv_dims(::HorizontalAverage, out_dims::Tuple) = tuple(last(out_dims))

macro horizontal_average(config_type, name)
    iex = generate_dv_interface(:HorizontalAverage, config_type, name)
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro horizontal_average(
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = generate_dv_interface(
        :HorizontalAverage,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro horizontal_average(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = quote
        $(generate_dv_interface(
            :HorizontalAverage,
            config_type,
            name,
            units,
            long_name,
            standard_name,
        ))
        $(generate_dv_function(
            :HorizontalAverage,
            config_type,
            [name],
            impl,
        ))
    end
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro horizontal_average_impl(
    impl,
    config_type,
    names...,
)
    iex = generate_dv_function(
        :HorizontalAverage,
        config_type,
        names,
        impl,
    )
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

"""
    ScalarDiagnostic

A reduction into a scalar value.
"""
abstract type ScalarDiagnostic <: DiagnosticVar end
dv_ScalarDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{ScalarDiagnostic}},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_dims(::ScalarDiagnostic, ::Any) = ()

macro scalar_diagnostic(config_type, name)
    iex = generate_dv_interface(:ScalarDiagnostic, config_type, name)
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro scalar_diagnostic(
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = generate_dv_interface(
        :ScalarDiagnostic,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro scalar_diagnostic(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = quote
        $(generate_dv_interface(
            :ScalarDiagnostic,
            config_type,
            name,
            units,
            long_name,
            standard_name,
        ))
        $(generate_dv_function(
            :ScalarDiagnostic,
            config_type,
            [name],
            impl,
        ))
    end
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

macro scalar_diagnostic_impl(
    impl,
    config_type,
    names...,
)
    iex = generate_dv_function(
        :ScalarDiagnostic,
        config_type,
        names,
        impl,
    )
    esc(MacroTools.prewalk(unblock, MacroTools.prewalk(rmlines, iex)))
end

include("atmos_diagnostic_funs.jl")