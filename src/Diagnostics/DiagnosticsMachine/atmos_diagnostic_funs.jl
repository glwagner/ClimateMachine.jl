using ..Atmos
using ..TurbulenceClosures
using ..TurbulenceConvection

# Method definitions for diagnostics collection for all the sub-models
# in `AtmosModel`.

dv_IntermediateValue(
    ::Type{AtmosConfigType},
    ::Union{Type{IntermediateValue}},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{AtmosConfigType},
    ::Union{Type{HorizontalAverage}},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{ScalarDiagnostic}},
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{AtmosConfigType},
    ::Union{Type{IntermediateValue}},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{AtmosConfigType},
    ::Union{Type{HorizontalAverage}},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{ScalarDiagnostic}},
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{AtmosConfigType},
    ::Union{Type{IntermediateValue}},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{AtmosConfigType},
    ::Union{Type{HorizontalAverage}},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{ScalarDiagnostic}},
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{AtmosConfigType},
    ::Union{Type{IntermediateValue}},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{AtmosConfigType},
    ::Union{Type{HorizontalAverage}},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{ScalarDiagnostic}},
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{AtmosConfigType},
    ::Union{Type{IntermediateValue}},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{AtmosConfigType},
    ::Union{Type{HorizontalAverage}},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{ScalarDiagnostic}},
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing

dv_IntermediateValue(
    ::Type{AtmosConfigType},
    ::Union{Type{IntermediateValue}},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_PointwiseDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{PointwiseDiagnostic}},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_HorizontalAverage(
    ::Type{AtmosConfigType},
    ::Union{Type{HorizontalAverage}},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing
dv_ScalarDiagnostic(
    ::Type{AtmosConfigType},
    ::Union{Type{ScalarDiagnostic}},
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::AbstractIntermediates,
) = nothing