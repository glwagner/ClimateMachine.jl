@pointwise_diagnostic(
    AtmosGCMConfigType,
    "u",
    "m s^-1",
    "zonal wind",
    "eastward_wind",
) do (atmos::AtmosModel, states, curr_time)
    u = states.prognostic.ρu[1] / states.prognostic.ρ
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "v",
    "m s^-1",
    "meridional wind",
    "northward_wind",
) do (atmos::AtmosModel, states, curr_time)
    v = states.prognostic.ρu[2] / states.prognostic.ρ
end


@pointwise_diagnostic(
    AtmosGCMConfigType,
    "w",
    "m s^-1",
    "vertical wind",
    "upward_air_velocity",
) do (atmos::AtmosModel, states, curr_time)
    w = states.prognostic.ρu[3] / states.prognostic.ρ
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "rho",
    "kg m^-3",
    "air density",
    "air_density",
) do (atmos::AtmosModel, states, curr_time)
    rho = states.prognostic.ρ
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "temp",
    "K",
    "air temperature",
    "air_temperature",
) do (atmos::AtmosModel, states, curr_time)
    temp = states.thermodynamic.temp
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "pres",
    "Pa",
    "air pressure",
    "air_pressure",
) do (atmos::AtmosModel, states, curr_time)
    pres = states.thermodynamic.pres
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "thd",
    "K",
    "dry potential temperature",
    "air_potential_temperature",
) do (atmos::AtmosModel, states, curr_time)
    thd = states.thermodynamic.θ_dry
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "et",
    "J kg^-1",
    "total specific energy",
    "specific_dry_energy_of_air",
) do (atmos::AtmosModel, states, curr_time)
    et = states.prognostic.ρe / states.prognostic.ρ
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "ei",
    "J kg^-1",
    "specific internal energy",
    "internal_energy",
) do (atmos::AtmosModel, states, curr_time)
    ei = states.thermodynamic.e_int
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "ht",
    "J kg^-1",
    "specific enthalpy based on total energy",
    "",
) do (atmos::AtmosModel, states, curr_time)
    ht = states.thermodynamic.h_tot
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "hi",
    "J kg^-1",
    "specific enthalpy based on internal energy",
    "atmosphere_enthalpy_content",
) do (atmos::AtmosModel, states, curr_time)
    hi = states.thermodynamic.h_int
end

#= TODO
@XXX_diagnostic(
    "vort",
    AtmosGCMConfigType,
    GridInterpolated,
    "s^-1",
    "vertical component of relative velocity",
    "atmosphere_relative_velocity",
) do (atmos::AtmosModel, states, curr_time)
end
=#

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "qt",
    "kg kg^-1",
    "mass fraction of total water in air (qv+ql+qi)",
    "mass_fraction_of_water_in_air",
) do (m::Union{EquilMoist, NonEquilMoist}, atmos::AtmosModel, states, curr_time)
    qt = states.prognostic.moisture.ρq_tot / states.prognostic.ρ
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "ql",
    "kg kg^-1",
    "mass fraction of liquid water in air",
    "mass_fraction_of_cloud_liquid_water_in_air",
) do (m::Union{EquilMoist, NonEquilMoist}, atmos::AtmosModel, states, curr_time)
    ql = states.thermodynamic.moisture.q_liq
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "qv",
    "kg kg^-1",
    "mass fraction of water vapor in air",
    "specific_humidity",
) do (m::Union{EquilMoist, NonEquilMoist}, atmos::AtmosModel, states, curr_time)
    qv = states.thermodynamic.moisture.q_vap
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "qi",
    "kg kg^-1",
    "mass fraction of ice in air",
    "mass_fraction_of_cloud_ice_in_air",
) do (m::Union{EquilMoist, NonEquilMoist}, atmos::AtmosModel, states, curr_time)
    qi = states.thermodynamic.moisture.q_ice
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "thv",
    "K",
    "virtual potential temperature",
    "virtual_potential_temperature",
) do (m::Union{EquilMoist, NonEquilMoist}, atmos::AtmosModel, states, curr_time)
    thv = states.thermodynamic.moisture.θ_vir
end

@pointwise_diagnostic(
    AtmosGCMConfigType,
    "thl",
    "K",
    "liquid-ice potential temperature",
    "",
) do (m::Union{EquilMoist, NonEquilMoist}, atmos::AtmosModel, states, curr_time)
    thl = states.thermodynamic.moisture.θ_liq_ice
end
