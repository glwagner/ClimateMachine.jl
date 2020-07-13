using ..Atmos
using ..TurbulenceClosures
using ..TurbulenceConvection

# Method definitions for diagnostics collection for all the sub-models
# in `AtmosModel`.

dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::MoistureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{PointwiseDiagnostic},
    ::MoistureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Type{HorizontalAverage},
    ::MoistureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{ScalarDiagnostic},
    ::MoistureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::PrecipitationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{PointwiseDiagnostic},
    ::PrecipitationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Type{HorizontalAverage},
    ::PrecipitationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{ScalarDiagnostic},
    ::PrecipitationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::RadiationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{PointwiseDiagnostic},
    ::RadiationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Type{HorizontalAverage},
    ::RadiationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{ScalarDiagnostic},
    ::RadiationModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::TracerModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{PointwiseDiagnostic},
    ::TracerModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Type{HorizontalAverage},
    ::TracerModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{ScalarDiagnostic},
    ::TracerModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::TurbulenceClosureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{PointwiseDiagnostic},
    ::TurbulenceClosureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Type{HorizontalAverage},
    ::TurbulenceClosureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{ScalarDiagnostic},
    ::TurbulenceClosureModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{ClimateMachineConfigType},
    ::Union{Type{IntermediateValue}},
    ::TurbulenceConvectionModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{PointwiseDiagnostic},
    ::TurbulenceConvectionModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{ClimateMachineConfigType},
    ::Type{HorizontalAverage},
    ::TurbulenceConvectionModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{ClimateMachineConfigType},
    ::Type{ScalarDiagnostic},
    ::TurbulenceConvectionModel,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing