#####
##### Tendency specification
#####

# Mass
ρ_tend(::AtmosModel, ::Flux1ˢᵗOrder) = (ρAdvect(),)

# Momentum
ρu_pressure(ref_state::HydrostaticState) = ρuHydroPress()
ρu_pressure(ref_state) = ρuNonHydroPress()

ρu_tend(m::AtmosModel, ::Flux1ˢᵗOrder) = (ρuAdvect(), ρu_pressure(m.ref_state))

# Energy
ρe_tend(::AtmosModel, ::Flux1ˢᵗOrder) = (ρeAdvect(), ρePressure())
