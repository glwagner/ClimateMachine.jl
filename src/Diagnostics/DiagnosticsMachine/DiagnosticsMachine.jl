"""
    DiagnosticsMachine

This module provides the infrastructure to extract diagnostics from a
ClimateMachine simulation. Two key abstractions are defined: diagnostic
variables and diagnostic groups. The `StdDiagnostics` module makes use of
these to define many standard variables and groups which may be used
directly by experiments. `DiagnosticsMachine` may be used by experiments
to define specialized variables and groups.
"""
module DiagnosticsMachine

export DiagnosticVar,
    dv_name,
    dv_attrib,
    dv_args,
    IntermediateValue,
    @intermediate_value,
    @intermediate_values,
    PointwiseDiagnostic,
    @pointwise_diagnostic,
    @pointwise_diagnostic_impl,
    HorizontalAverage,
    @horizontal_average,
    @horizontal_average_impl,
    ScalarDiagnostic,
    @scalar_diagnostic,
    @scalar_diagnostic_impl,
    States,
    DiagnosticsGroup,
    @diagnostics_group,
    DiagnosticsGroupParams

using CUDA
using Dates
using InteractiveUtils
using KernelAbstractions
using MacroTools
using MPI
using OrderedCollections
using Printf

using ..Atmos
using ..BalanceLaws
using ..ConfigTypes
using ..DGMethods
using ..GenericCallbacks
using ..Mesh.Interpolation
using ..MPIStateArrays
using ..Spectra
using ..TicToc
using ..VariableTemplates
using ..Writers

using CLIMAParameters
using CLIMAParameters.Planet: planet_radius

# Container to store simulation information necessary for all
# diagnostics groups.
Base.@kwdef mutable struct Diagnostic_Settings
    mpicomm::MPI.Comm = MPI.COMM_WORLD
    param_set::Union{Nothing, AbstractParameterSet} = nothing
    dg::Union{Nothing, DGModel} = nothing
    Q::Union{Nothing, MPIStateArray} = nothing
    starttime::Union{Nothing, String} = nothing
    output_dir::Union{Nothing, String} = nothing
end
const Settings = Diagnostic_Settings()

include("variables.jl")
include("groups.jl")

const AllDiagnosticVars =
    OrderedDict{
        Type{<:ClimateMachineConfigType},
        OrderedDict{String, Type{<:DiagnosticVar}},
    }()
AllDiagnosticVars[ClimateMachineConfigType] =
    OrderedDict{String, Type{<:DiagnosticVar}}()
function add_all_dvar_dicts(T::DataType)
    for t in subtypes(T)
        AllDiagnosticVars[t] =
            OrderedDict{String, Type{<:DiagnosticVar}}()
        add_all_dvar_dicts(t)
    end
end
add_all_dvar_dicts(ClimateMachineConfigType)
#dvar = AllDiagnosticVars[getfield(ConfigTypes, config_type)][name]

"""
    init(mpicomm, param_set, dg, Q, starttime, output_dir)

Save the parameters into `Settings`, a container for simulation
information necessary for all diagnostics groups.
"""
function init(
    mpicomm::MPI.Comm,
    param_set::AbstractParameterSet,
    dg::DGModel,
    Q::MPIStateArray,
    starttime::String,
    output_dir::String,
)
    Settings.mpicomm = mpicomm
    Settings.param_set = param_set
    Settings.dg = dg
    Settings.Q = Q
    Settings.starttime = starttime
    Settings.output_dir = output_dir
end

end # module DiagnosticsMachine
