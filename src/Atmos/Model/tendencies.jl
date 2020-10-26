#####
##### Tendency specification
#####

# Mass
ρ_sources(s) = ()
ρ_tend(::AtmosModel, ::Flux1ˢᵗOrder) = (ρAdvect(),)
ρ_tend(m::AtmosModel, ::NonConservative) = (ρ_sources.(m.source...)...,)

# Momentum
ρu_sources(s) = ()
ρu_pressure(ref_state::HydrostaticState) = ρuHydroPress()
ρu_pressure(ref_state) = ρuNonHydroPress()
ρu_tend(m::AtmosModel, ::Flux1ˢᵗOrder) = (ρuAdvect(), ρu_pressure(m.ref_state))
ρu_tend(m::AtmosModel, ::NonConservative) = (ρu_sources.(m.source...)...,)


# Energy
ρe_sources(s) = ()
ρe_sources(s::Subsidence) = ρeSubsidence(s)
ρe_tend(::AtmosModel, ::Flux1ˢᵗOrder) = (ρeAdvect(), ρePressure())
ρe_tend(m::AtmosModel, ::NonConservative) = (ρe_sources.(m.source...)...,)

# Moisture
ρq_tot_sources(s) = ()
ρq_tot_sources(s::Subsidence) = ρq_totSubsidence(s)
ρq_tot_tend(m::AtmosModel, ::NonConservative) = (ρq_tot_sources.(m.source...)...,)
