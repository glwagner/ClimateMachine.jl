##### Tendencies

#####
##### First order fluxes
#####

function Σfluxes(fluxes, m, state, aux, t, ts, direction)
    return sum(ntuple(length(fluxes)) do i
        _flux_(fluxes[i], m, state, aux, t, ts, direction)
    end)
end

#####
##### Second order fluxes
#####

function Σfluxes(fluxes, m, state, aux, t, ts, diffusive, hyperdiffusive)
    return sum(ntuple(length(fluxes)) do i
        _flux_(fluxes[i], m, state, aux, t, ts, diffusive, hyperdiffusive)
    end)
end

#####
##### Non-conservative sources
#####

function Σsources(sources, m, state, aux, t, ts, direction, diffusive)
    return sum(ntuple(length(sources)) do i
        _source_(sources[i], m, state, aux, t, ts, direction, diffusive)
    end)
end

#####
##### Tendency specification
#####

# Mass
ρ_moisture(::MoistModel) = (ρMoisture(),)
ρ_moisture(::MoistureModel) = ()

ρ_tend(::AtmosModel, ::Flux1ˢᵗOrder) = (ρAdvect(),)

ρ_tend(m::AtmosModel, ::Flux2ⁿᵈOrder) = (ρ_moisture(m)...,)

ρ_tend(m::AtmosModel, ::NonConservative) = ()

# Momentum
ρu_pressure(ref_state::HydrostaticState) = ρuHydroPress()
ρu_pressure(ref_state) = ρuNonHydroPress()

ρu_tend(m::AtmosModel, ::Flux1ˢᵗOrder) = (ρuAdvect(), ρu_pressure(m.ref_state))

ρu_tend(m::AtmosModel, ::Flux2ⁿᵈOrder) = (ρuShear(),)

ρu_tend(m::AtmosModel, ::NonConservative) = ()

# Energy
ρe_radiation(::DYCOMSRadiation) = (ρeDycomsRadiation(),)
ρe_radiation(::RadiationModel) = ()

ρe_tend(m::AtmosModel, ::Flux1ˢᵗOrder) =
    (ρeAdvect(), ρePressure(), ρe_radiation(m.radiation)...,)

ρe_tend(m::AtmosModel, ::Flux2ⁿᵈOrder) =
    (ρeShear(), ρeProduction(), subsidence(m)...)

ρe_tend(m::AtmosModel, ::NonConservative) = ()

# Total specific humidity
ρq_tot_tend(m::AtmosModel, ::AbstractSource) = ρq_tot_tend(m.moisture)
ρq_tot_tend(m::DryModel, ::Flux1ˢᵗOrder) = ()
ρq_tot_tend(m::MoistModel, ::Flux1ˢᵗOrder) = (ρq_totAdvect(),)

ρq_tot_tend(m::DryModel, ::Flux2ⁿᵈOrder) = ()
ρq_tot_tend(m::MoistModel, ::Flux2ⁿᵈOrder) = (ρq_totDiffusion(),)

ρq_tot_tend(m::MoistureModel, ::NonConservative) = (ρq_tot_subsidence(m)...,)

# Liquid specific humidity
ρq_liq_tend(m::AtmosModel, ::AbstractSource) = ρq_liq_tend(m.moisture)
ρq_liq_tend(m::MoistureModel, ::Flux1ˢᵗOrder) = ()
ρq_liq_tend(m::NonEquilMoist, ::Flux1ˢᵗOrder) = (ρq_liqAdvect(),)

ρq_liq_tend(m::MoistureModel, ::Flux2ⁿᵈOrder) = ()
ρq_liq_tend(m::NonEquilMoist, ::Flux2ⁿᵈOrder) = (ρq_liqDiffusion(),)

ρq_liq_tend(m::MoistureModel, ::NonConservative) = (ρq_liqCreateClouds(),)

# Ice specific humidity
ρq_ice_tend(m::AtmosModel, ::AbstractSource) = ρq_ice_tend(m.moisture)
ρq_ice_tend(m::MoistureModel, ::Flux1ˢᵗOrder) = ()
ρq_ice_tend(m::NonEquilMoist, ::Flux1ˢᵗOrder) = (ρq_iceAdvect(),)

ρq_ice_tend(m::MoistureModel, ::Flux2ⁿᵈOrder) = ()
ρq_ice_tend(m::NonEquilMoist, ::Flux2ⁿᵈOrder) = (ρq_iceDiffusion(),)

ρq_ice_tend(m::MoistureModel, ::NonConservative) = ()

# Interaction/multi-physics terms
# Example use:
subsidence(m::AtmosModel) = Subsidence{eltype(m)}(2.0)
""" Apply subsidence to ρe and ρq_tot equations """
subsidence(m::AtmosModel) = ()
ρe_subsidence(m::AtmosModel) = ρeSubsidence(subsidence(m))
ρq_tot_subsidence(m::AtmosModel) = ρq_totSubsidence(subsidence(m))

""" Create clouds in ρq_liq and ρq_ice equations """
CreateClouds(m::AtmosModel) = ()
ρq_liqCreateClouds(m::AtmosModel) = ρq_liqCreateClouds(CreateClouds(m))
ρq_iceCreateClouds(m::AtmosModel) = ρq_iceCreateClouds(CreateClouds(m))

