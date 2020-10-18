##### Moisture

#####
##### First order fluxes
#####

struct ρq_totAdvect <: Flux1ˢᵗOrder end
function _flux_(::ρq_totAdvect, m::AtmosModel, state, aux, t, ts, direction)
    return (state.ρu / state.ρ) * state.moisture.ρq_tot
end

struct ρq_liqAdvect <: Flux1ˢᵗOrder end
function _flux_(::ρq_liqAdvect, m::AtmosModel, state, aux, t, ts, direction)
    return (state.ρu / state.ρ) * state.moisture.ρq_liq
end

struct ρq_iceAdvect <: Flux1ˢᵗOrder end
function _flux_(::ρq_iceAdvect, m::AtmosModel, state, aux, t, ts, direction)
    return (state.ρu / state.ρ) * state.moisture.ρq_ice
end

#####
##### Second order fluxes
#####

struct ρq_totDiffusion <: Flux2ⁿᵈOrder end # name?
function _flux_(::ρq_totDiffusion, m::AtmosModel, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(m, state, diffusive, aux, t)
    d_q_tot = (-D_t) .* diffusive.moisture.∇q_tot
    return d_q_tot * state.ρ
end

struct ρq_liqDiffusion <: Flux2ⁿᵈOrder end # name?
function _flux_(::ρq_liqDiffusion, m::AtmosModel, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(m, state, diffusive, aux, t)
    d_q_liq = (-D_t) .* diffusive.moisture.∇q_liq
    return d_q_liq * state.ρ
end

struct ρq_iceDiffusion <: Flux2ⁿᵈOrder end # name?
function _flux_(::ρq_iceDiffusion, m::AtmosModel, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(m, state, diffusive, aux, t)
    d_q_ice = (-D_t) .* diffusive.moisture.∇q_ice
    return d_q_ice * state.ρ
end

#####
##### Sources
#####

const ρq_totSubsidence = Subsidence
function _source_(s::ρq_totSubsidence, atmos, state, aux, t, ts, direction, diffusive)
    z = altitude(atmos, aux)
    w_sub = subsidence_velocity(s, z)
    k̂ = vertical_unit_vector(atmos, aux)
    return - state.ρ * w_sub * dot(k̂, diffusive.moisture.∇q_tot)
end
ρq_totSubsidence(s::Subsidence{FT}) where {FT} = ρq_totSubsidence{FT}(s.D)
ρq_totSubsidence(::Tuple{}) = ()

struct CreateClouds <: Source end

const ρq_liqCreateClouds = CreateClouds
function _source_(s::ρq_liqCreateClouds, atmos, state, aux, t, ts, direction, diffusive)
    # get current temperature and phase partition
    ts = recover_thermo_state(atmos, state, aux)
    q = PhasePartition(ts)
    T = air_temperature(ts)

    # phase partition corresponding to the current T and q.tot
    # (this is not the same as phase partition from saturation adjustment)
    ts_eq = PhaseEquil_ρTq(atmos.param_set, state.ρ, T, q.tot)
    q_eq = PhasePartition(ts_eq)

    # cloud condensate as relaxation source terms
    S_q_liq = conv_q_vap_to_q_liq_ice(atmos.param_set.microphys.liq, q_eq, q)
    return state.ρ * S_q_liq
end
ρq_liqCreateClouds(::CreateClouds) = (ρq_liqCreateClouds(),)
ρq_liqCreateClouds(::Tuple{}) = ()

const ρq_iceCreateClouds = CreateClouds
function _source_(s::ρq_iceCreateClouds, atmos, state, aux, t, ts, direction, diffusive)
    # get current temperature and phase partition
    ts = recover_thermo_state(atmos, state, aux)
    q = PhasePartition(ts)
    T = air_temperature(ts)

    # phase partition corresponding to the current T and q.tot
    # (this is not the same as phase partition from saturation adjustment)
    ts_eq = PhaseEquil_ρTq(atmos.param_set, state.ρ, T, q.tot)
    q_eq = PhasePartition(ts_eq)

    # cloud condensate as relaxation source terms
    S_q_ice = conv_q_vap_to_q_liq_ice(atmos.param_set.microphys.ice, q_eq, q)
    return state.ρ * S_q_ice
end
ρq_iceCreateClouds(::CreateClouds) = (ρq_iceCreateClouds(),)
ρq_iceCreateClouds(::Tuple{}) = ()

