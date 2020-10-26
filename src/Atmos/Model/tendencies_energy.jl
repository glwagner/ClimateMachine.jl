##### Energy tendencies

#####
##### First order fluxes
#####

struct ρeAdvect <: Flux1ˢᵗOrderTendency end
function _flux_(::ρeAdvect, m, state, aux, t, ts, direction)
    return (state.ρu / state.ρ) * state.ρe
end

struct ρePressure <: Flux1ˢᵗOrderTendency end
function _flux_(::ρePressure, m, state, aux, t, ts, direction)
    return state.ρu / state.ρ * air_pressure(ts)
end

#####
##### Non-conservative sources
#####

const ρeSubsidence = Subsidence
function _source_(s::ρeSubsidence, m, state, aux, t, ts, direction, diffusive)
    z = altitude(m, aux)
    w_sub = subsidence_velocity(s, z)
    k̂ = vertical_unit_vector(m, aux)
    return -state.ρ * w_sub * dot(k̂, diffusive.∇h_tot)
end
ρeSubsidence(s::Subsidence{FT}) where {FT} = ρeSubsidence{FT}(s.D)
