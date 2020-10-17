##### Energy tendencies

struct ρeAdvect <: Flux1ˢᵗOrderTendency end
function _flux_(::ρeAdvect, m, state, aux, t, ts, direction)
    return (state.ρu / state.ρ) * state.ρe
end

struct ρePressure <: Flux1ˢᵗOrderTendency end
function _flux_(::ρePressure, m, state, aux, t, ts, direction)
    return state.ρu / state.ρ * air_pressure(ts)
end
