##### Multi-physics tendencies

struct Subsidence{FT} <: NonConservativeTendency
    D::FT
end
subsidence_velocity(subsidence::Subsidence{FT}, z::FT) where {FT} =
    -subsidence.D * z
