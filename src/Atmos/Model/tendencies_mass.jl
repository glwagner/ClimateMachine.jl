##### Mass tendencies

struct ρAdvect <: Flux1ˢᵗOrderTendency end
function _flux_(::ρAdvect, m, state, aux, t, ts, direction)
    return state.ρu
end
