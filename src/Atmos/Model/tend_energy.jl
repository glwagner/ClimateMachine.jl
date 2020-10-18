##### Energy

#####
##### First order fluxes
#####

struct ρeAdvect <: Flux1ˢᵗOrder end
function _flux_(::ρeAdvect, m, state, aux, t, ts, direction)
    return (state.ρu / state.ρ) * state.ρe
end

struct ρePressure <: Flux1ˢᵗOrder end
function _flux_(::ρePressure, m, state, aux, t, ts, direction)
    return (state.ρu / state.ρ) * air_pressure(ts)
end

struct ρeDycomsRadiation <: Flux1ˢᵗOrder end
function _flux_(::ρeDycomsRadiation, m, state, aux, t, ts, direction)
    FT = eltype(flux)
    z = altitude(atmos, aux)
    Δz_i = max(z - m.z_i, -zero(FT))
    # Constants
    upward_flux_from_cloud = m.F_0 * exp(-aux.∫dnz.radiation.attenuation_coeff)
    upward_flux_from_sfc = m.F_1 * exp(-aux.∫dz.radiation.attenuation_coeff)
    free_troposphere_flux =
        m.ρ_i *
        FT(cp_d(atmos.param_set)) *
        m.D_subsidence *
        m.α_z *
        cbrt(Δz_i) *
        (Δz_i / 4 + m.z_i)
    F_rad =
        upward_flux_from_sfc + upward_flux_from_cloud + free_troposphere_flux
    ẑ = vertical_unit_vector(atmos, aux)
    return F_rad * ẑ
end

#####
##### Second order fluxes
#####

struct ρeShear <: Flux2ⁿᵈOrder end # name?
function _flux_(::ρeShear, m, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(m, state, diffusive, aux, t)
    return τ * state.ρu
end

struct ρeProduction <: Flux2ⁿᵈOrder end # name?
function _flux_(::ρeProduction, m, state, aux, t, ts, diffusive, hyperdiffusive)
    ν, D_t, τ = turbulence_tensors(m, state, diffusive, aux, t)
    d_h_tot = -D_t .* diffusive.∇h_tot
    return d_h_tot * state.ρ
end

#####
##### Sources
#####

const ρeSubsidence = Subsidence
function _source_(s::ρeSubsidence, atmos, state, aux, t, ts, direction, diffusive)
    z = altitude(atmos, aux)
    w_sub = subsidence_velocity(s, z)
    k̂ = vertical_unit_vector(atmos, aux)
    return - state.ρ * w_sub * dot(k̂, diffusive.∇h_tot)
end
ρeSubsidence(s::Subsidence{FT}) where {FT} = ρeSubsidence{FT}(s.D)
ρeSubsidence(::Tuple{}) = ()
