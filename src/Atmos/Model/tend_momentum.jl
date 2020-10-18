##### Momentum

#####
##### First order fluxes
#####

struct ρuAdvect <: Flux1ˢᵗOrder end
function _flux_(::ρuAdvect, atmos, state, aux, t, ts, direction)
    return state.ρu .* (state.ρu / state.ρ)'
end

struct ρuHydroPress <: Flux1ˢᵗOrder end
function _flux_(::ρuHydroPress, atmos, state, aux, t, ts, direction)
    return (air_pressure(ts) - aux.ref_state.p) * I
end

struct ρuNonHydroPress <: Flux1ˢᵗOrder end
function _flux_(::ρuNonHydroPress, atmos, state, aux, t, ts, direction)
    return air_pressure(ts) * I
end

#####
##### Second order fluxes
#####

struct ρuShear <: Flux2ⁿᵈOrder end
function _flux_(::ρuShear, atmos, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(atmos, state, diffusive, aux, t)
    return τ * state.ρ
end

struct ρuMoisture <: Flux2ⁿᵈOrder end # name?
function _flux_(::ρuMoisture, atmos, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(atmos, state, diffusive, aux, t)
    d_q_tot = (-D_t) .* diffusive.moisture.∇q_tot
    return d_q_tot .* state.ρu'
end

#####
##### Sources
#####

struct ρuGravity <: NonConservative end
function _source_(::ρuGravity, atmos, state, aux, t, ts, direction, diffusive)
    if m.ref_state isa HydrostaticState
        return -(state.ρ - aux.ref_state.ρ) * aux.orientation.∇Φ
    else
        return -state.ρ * aux.orientation.∇Φ
    end
end

struct ρuCoriolis <: NonConservative end
function _source_(::ρuCoriolis, atmos, state, aux, t, ts, direction, diffusive)
    FT = eltype(state)
    _Omega::FT = Omega(m.param_set)
    # note: this assumes a SphericalOrientation
    return - SVector(0, 0, 2 * _Omega) × state.ρu
end

struct GeostrophicForcing{FT} <: NonConservative
    f_coriolis::FT
    u_geostrophic::FT
    v_geostrophic::FT
end
function _source_(::ρuGeostrophicForcing, atmos, state, aux, t, ts, direction, diffusive)
    u_geo = SVector(s.u_geostrophic, s.v_geostrophic, 0)
    ẑ = vertical_unit_vector(atmos, aux)
    fkvector = s.f_coriolis * ẑ
    return - fkvector × (state.ρu .- state.ρ * u_geo)
end

"""
    ρuRayleighSponge{FT} <: Source

Rayleigh Damping (Linear Relaxation) for top wall momentum components
Assumes laterally periodic boundary conditions for LES flows. Momentum components
are relaxed to reference values (zero velocities) at the top boundary.
"""
struct ρuRayleighSponge{FT} <: NonConservative
    "Maximum domain altitude (m)"
    z_max::FT
    "Altitude at with sponge starts (m)"
    z_sponge::FT
    "Sponge Strength 0 ⩽ α_max ⩽ 1"
    α_max::FT
    "Relaxation velocity components"
    u_relaxation::SVector{3, FT}
    "Sponge exponent"
    γ::FT
end
function _source_(s::ρuRayleighSponge, atmos, state, aux, t, ts, direction, diffusive)
    z = altitude(atmos, aux)
    if z >= s.z_sponge
        FT = eltype(state)
        r = (z - s.z_sponge) / (s.z_max - s.z_sponge)
        β_sponge = s.α_max * sinpi(r / 2)^s.γ
        return - β_sponge * (state.ρu .- state.ρ * s.u_relaxation)
    else
        return FT(0)
    end
end
