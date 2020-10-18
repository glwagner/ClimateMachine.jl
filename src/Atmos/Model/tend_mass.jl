##### Mass

#####
##### First order fluxes
#####

struct ρAdvect <: Flux1ˢᵗOrder end
function _flux_(::ρAdvect, m, state, aux, t, ts, direction)
    return state.ρu
end

#####
##### Second order fluxes
#####

struct ρMoisture <: Flux2ⁿᵈOrder end # name?
function _flux_(::ρMoisture, m, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(m, state, diffusive, aux, t)
    d_q_tot = (-D_t) .* diffusive.moisture.∇q_tot
    return d_q_tot * state.ρ
end

