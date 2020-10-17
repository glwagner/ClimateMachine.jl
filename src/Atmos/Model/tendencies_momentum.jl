##### Momentum tendencies

struct ρuAdvect <: Flux1ˢᵗOrderTendency end
function _flux_(::ρuAdvect, m, state, aux, t, ts, direction)
    return state.ρu .* (state.ρu / state.ρ)'
end

struct ρuHydroPress <: Flux1ˢᵗOrderTendency end
function _flux_(::ρuHydroPress, m, state, aux, t, ts, direction)
    return (air_pressure(ts) - aux.ref_state.p) * I
end

struct ρuNonHydroPress <: Flux1ˢᵗOrderTendency end
function _flux_(::ρuNonHydroPress, m, state, aux, t, ts, direction)
    return air_pressure(ts) * I
end
