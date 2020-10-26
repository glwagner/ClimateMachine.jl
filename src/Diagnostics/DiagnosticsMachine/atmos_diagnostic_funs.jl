using ..Atmos
using ..TurbulenceClosures
using ..TurbulenceConvection

# Method definitions for diagnostics collection for all the sub-models
# in `AtmosModel`.

dv_IntermediateValue(
    ::AtmosConfigType,
    ::Union{IntermediateValue},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::Union{PointwiseDiagnostic},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::Union{HorizontalAverage},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::Union{ScalarDiagnostic},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::AtmosConfigType,
    ::Union{IntermediateValue},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::Union{PointwiseDiagnostic},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::Union{HorizontalAverage},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::Union{ScalarDiagnostic},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::AtmosConfigType,
    ::Union{IntermediateValue},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::Union{PointwiseDiagnostic},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::Union{HorizontalAverage},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::Union{ScalarDiagnostic},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::AtmosConfigType,
    ::Union{IntermediateValue},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::Union{PointwiseDiagnostic},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::Union{HorizontalAverage},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::Union{ScalarDiagnostic},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::AtmosConfigType,
    ::Union{IntermediateValue},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::Union{PointwiseDiagnostic},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::Union{HorizontalAverage},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::Union{ScalarDiagnostic},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::AtmosConfigType,
    ::Union{IntermediateValue},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::Union{PointwiseDiagnostic},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::Union{HorizontalAverage},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::Union{ScalarDiagnostic},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing