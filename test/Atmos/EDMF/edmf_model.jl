#### EDMF model

Base.@kwdef struct Updraft{FT <: AbstractFloat} <: BalanceLaw end

Base.@kwdef struct EDMF{FT <: AbstractFloat, N, UP, TCD} <: TurbulenceConvectionModel
    "Updrafts"
    updraft::UP
    "Turbconv data"
    turbconv_data::TCD
end

function EDMF(
    FT,
    N_up;
    updraft = ntuple(i -> Updraft{FT}(), N_up),
    turbconv_data=turbconv_data,
)
    args = (updraft,turbconv_data)
    return EDMF{FT, N_up, typeof.(args)...}(args...)
end

import ClimateMachine.TurbulenceConvection: turbconv_sources, turbconv_bcs

struct EDMFBCs <: TurbConvBC end
n_updrafts(m::EDMF{FT, N_up}) where {FT, N_up} = N_up
turbconv_filters(m::EDMF) = ("turbconv.updraft",)
turbconv_sources(m::EDMF) = (turbconv_source!,)
turbconv_bcs(::EDMF) = EDMFBCs()
