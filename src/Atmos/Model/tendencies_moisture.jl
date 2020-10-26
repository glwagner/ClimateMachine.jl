##### Moisture tendencies

#####
##### Non-conservative sources
#####

const ρq_totSubsidence = Subsidence
function _source_(
    s::ρq_totSubsidence,
    m,
    state,
    aux,
    t,
    ts,
    direction,
    diffusive,
)
    z = altitude(m, aux)
    w_sub = subsidence_velocity(s, z)
    k̂ = vertical_unit_vector(m, aux)
    return -state.ρ * w_sub * dot(k̂, diffusive.moisture.∇q_tot)
end
ρq_totSubsidence(s::Subsidence{FT}) where {FT} = ρq_totSubsidence{FT}(s.D)
