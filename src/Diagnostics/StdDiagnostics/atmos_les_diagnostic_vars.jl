@intermediate_value(AtmosLESConfigType, "D_t") do (
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    _, D_t, _ = turbulence_tensors(
        atmos,
        states.prognostic,
        states.gradient_flux,
        states.auxiliary,
        curr_time,
    )
end

@horizontal_average(AtmosLESConfigType, "u", "m s^-1", "x-velocity", "")
@horizontal_average(AtmosLESConfigType, "v", "m s^-1", "y-velocity", "")
@horizontal_average(AtmosLESConfigType, "w", "m s^-1", "z-velocity", "")

@horizontal_average(AtmosLESConfigType, "u_u")
@horizontal_average(AtmosLESConfigType, "v_v")
@horizontal_average(AtmosLESConfigType, "w_w")
@horizontal_average(AtmosLESConfigType, "w_w_w")
@horizontal_average(AtmosLESConfigType, "w_u")
@horizontal_average(AtmosLESConfigType, "w_v")

@horizontal_average_impl(
    AtmosLESConfigType,
    "u",
    "v",
    "w",
    "u_u",
    "v_v",
    "w_w",
    "w_w_w",
    "w_u",
    "w_v",
) do (atmos::AtmosModel, states::States, curr_time, intermediates)
    u = states.prognostic.ρu[1]
    v = states.prognostic.ρu[2]
    w = states.prognostic.ρu[3]
    u_u = u^2 / states.prognostic.ρ
    v_v = v^2 / states.prognostic.ρ
    w_w = w^2 / states.prognostic.ρ
    w_w_w = w^3 / states.prognostic.ρ^2
    w_u = w * u / states.prognostic.ρ
    w_v = w * v / states.prognostic.ρ
end

@horizontal_average(AtmosLESConfigType, "avg_rho", "kg m^-3", "air density", "air_density")
@horizontal_average(AtmosLESConfigType, "rho", "kg m^-3", "air density", "air_density")

@horizontal_average(AtmosLESConfigType, "w_rho")

@horizontal_average_impl(
    AtmosLESConfigType,
    "avg_rho",
    "rho",
    "w_rho",
) do (
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    avg_rho = states.prognostic.ρ
    rho = states.prognostic.ρ * states.prognostic.ρ
    wrho = states.prognostic.ρu[3] * states.prognostic.ρ
end

@horizontal_average(AtmosLESConfigType, "temp", "K", "air temperature", "air_temperature")
@horizontal_average(AtmosLESConfigType, "pres", "Pa", "air pressure", "air_pressure")
@horizontal_average(AtmosLESConfigType, "thd", "K", "dry potential temperature", "air_potential_temperature")
@horizontal_average(AtmosLESConfigType, "et", "J kg^-1", "total specific energy", "specific_dry_energy_of_air")
@horizontal_average(AtmosLESConfigType, "ei", "J kg^-1", "specific internal energy", "internal_energy")
@horizontal_average(AtmosLESConfigType, "ht", "J kg^-1", "specific enthalpy based on total energy", "")
@horizontal_average(AtmosLESConfigType, "hi", "J kg^-1", "specific enthalpy based on internal energy", "atmosphere_enthalpy_content")
@horizontal_average(AtmosLESConfigType, "w_ht_sgs", "kg kg^-1 m s^-1", "vertical sgs flux of total specific enthalpy", "")

@horizontal_average(AtmosLESConfigType, "ei_ei")
@horizontal_average(AtmosLESConfigType, "w_thd")
@horizontal_average(AtmosLESConfigType, "w_ei")

@horizontal_average_impl(
    AtmosLESConfigType,
    "temp",
    "pres",
    "et",
    "ei",
    "ht",
    "hi",
    "w_ht_sgs",
    "ei_ei",
    "w_thd",
    "w_ei",
) do (
    atmos::AtmosModel,
    states::States,
    curr_time,
    intermediates,
)
    temp = intermediates.temp * states.prognostic.ρ
    pres = intermediates.pres * states.prognostic.ρ
    thd = intermediates.θ_dry * states.prognostic.ρ
    et = states.prognostic.ρe
    ei = intermediates.e_int * states.prognostic.ρ
    ht = intermediates.h_tot * states.prognostic.ρ
    hi = intermediates.h_int * states.prognostic.ρ
    d_h_tot = -(intermediates.D_t) .* states.gradient_flux.∇h_tot
    w_ht_sgs = d_h_tot[end] * state.prognostic.ρ
    eiei = ei * intermediates.e_int
    wthd = states.prognostic.ρu[3] * intermediates.θ_dry
    wei = states.prognostic.ρu[3] * intermediates.e_int
end

@horizontal_average(AtmosLESConfigType, "qt", "kg kg^-1", "mass fraction of total water in air (qv+ql+qi)", "mass_fraction_of_water_in_air")
@horizontal_average(AtmosLESConfigType, "ql", "kg kg^-1", "mass fraction of liquid water in air", "mass_fraction_of_cloud_liquid_water_in_air")
@horizontal_average(AtmosLESConfigType, "qi", "kg kg^-1", "mass fraction of ice in air", "mass_fraction_of_cloud_ice_in_air")
@horizontal_average(AtmosLESConfigType, "qv", "kg kg^-1", "mass fraction of water vapor in air", "specific_humidity")
@horizontal_average(AtmosLESConfigType, "thv", "K", "virtual potential temperature", "virtual_potential_temperature")
@horizontal_average(AtmosLESConfigType, "thl", "K", "liquid-ice potential temperature", "")
@horizontal_average(AtmosLESConfigType, "w_qt_sgs", "kg kg^-1 m s^-1", "vertical sgs flux of total specific humidity", "")

@horizontal_average(AtmosLESConfigType, "qt_qt")
@horizontal_average(AtmosLESConfigType, "thl_thl")
@horizontal_average(AtmosLESConfigType, "w_qt")
@horizontal_average(AtmosLESConfigType, "w_ql")
@horizontal_average(AtmosLESConfigType, "w_qi")
@horizontal_average(AtmosLESConfigType, "w_qv")
@horizontal_average(AtmosLESConfigType, "w_thv")
@horizontal_average(AtmosLESConfigType, "w_thl")
@horizontal_average(AtmosLESConfigType, "qt_thl")
@horizontal_average(AtmosLESConfigType, "qt_ei")

@horizontal_average_impl(
    AtmosLESConfigType,
    "qt",
    "ql",
    "qi",
    "qv",
    "thv",
    "thl",
    "w_qt_sgs",
    "qt_qt",
    "thl_thl",
    "w_qt",
    "w_ql",
    "w_qi",
    "w_qv",
    "w_thv",
    "w_thl",
    "qt_thl",
    "qt_ei",
) do (
    m::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    intermediates
)
    qt = states.prognostic.moisture.ρq_tot
    ql = intermediates.moisture.q_liq * states.prognostic.ρ # TODO: `moisture.`?
    qi = intermediates.moisture.q_ice * states.prognostic.ρ
    qv = intermediates.moisture.q_vap * states.prognostic.ρ
    thv = intermediates.moisture.θ_vir * states.prognostic.ρ
    thl = intermediates.moisture.θ_liq_ice * states.prognostic.ρ
    d_q_tot = (-intermediates.D_t) .* states.gradient_flux.moisture.∇q_tot
    w_qt_sgs = d_q_tot[end] * states.prognostic.ρ
    qt_qt = qt * (states.prognostic.moisture.ρq_tot / states.prognostic.ρ)
    thl_thl = thl * intermediates.moisture.θ_liq_ice
    w = states.prognostic.ρu[3] / states.prognostic.ρ
    w_qt = states.prognostic.ρu[3] * qt / states.prognostic.ρ
    w_ql = states.prognostic.ρu[3] * intermediates.moisture.q_liq
    w_qi = states.prognostic.ρu[3] * intermediates.moisture.q_ice
    w_qv = states.prognostic.ρu[3] * intermediates.moisture.q_vap
    w_thv = states.prognostic.ρu[3] * intermediates.moisture.θ_vir
    w_thl = states.prognostic.ρu[3] * intermediates.moisture.θ_liq_ice
    qt_thl = qt * intermediates.moisture.θ_liq_ice
    qt_ei = qt * intermediates.e_int
end

#= TODO
    Variables["cld_frac"] = DiagnosticVariable(
        "cld_frac",
        diagnostic_var_attrib(
            "",
            "cloud fraction",
            "cloud_area_fraction_in_atmosphere_layer",
        ),
    )
    Variables["cld_cover"] = DiagnosticVariable(
        "cld_cover",
        diagnostic_var_attrib("", "cloud cover", "cloud_area_fraction"),
    )
    Variables["cld_top"] = DiagnosticVariable(
        "cld_top",
        diagnostic_var_attrib("m", "cloud top", "cloud_top_altitude"),
    )
    Variables["cld_base"] = DiagnosticVariable(
        "cld_base",
        diagnostic_var_attrib("m", "cloud base", "cloud_base_altitude"),
    )
    Variables["lwp"] = DiagnosticVariable(
        "lwp",
        diagnostic_var_attrib(
            "kg m^-2",
            "liquid water path",
            "atmosphere_mass_content_of_cloud_condensed_water",
        ),
    )
    Variables["core_frac"] = DiagnosticVariable(
        "core_frac",
        diagnostic_var_attrib("", "cloud core fraction", ""),
    )
    Variables["u_core"] = DiagnosticVariable(
        "u_core",
        diagnostic_var_attrib("m s^-1", "cloud core x-velocity", ""),
    )
    Variables["v_core"] = DiagnosticVariable(
        "v_core",
        diagnostic_var_attrib("m s^-1", "cloud core y-velocity", ""),
    )
    Variables["w_core"] = DiagnosticVariable(
        "w_core",
        diagnostic_var_attrib("m s^-1", "cloud core z-velocity", ""),
    )
    Variables["avg_rho_core"] = DiagnosticVariable(
        "avg_rho_core",
        diagnostic_var_attrib("kg m^-3", "cloud core air density", ""),
    )
    Variables["rho_core"] = DiagnosticVariable(
        "rho_core",
        diagnostic_var_attrib("kg m^-3", "cloud core (density-averaged) air density", ""),
    )
    Variables["qt_core"] = DiagnosticVariable(
        "qt_core",
        diagnostic_var_attrib("kg m^-3", "cloud core total specific humidity", ""),
    )
    Variables["ql_core"] = DiagnosticVariable(
        "ql_core",
        diagnostic_var_attrib("kg m^-3", "cloud core liquid water specific humidity", ""),
    )
    Variables["thv_core"] = DiagnosticVariable(
        "thv_core",
        diagnostic_var_attrib("K", "cloud core virtual potential temperature", ""),
    )
    Variables["thl_core"] = DiagnosticVariable(
        "thl_core",
        diagnostic_var_attrib("K", "cloud core liquid-ice potential temperature", ""),
    )
    Variables["ei_core"] = DiagnosticVariable(
        "ei_core",
        diagnostic_var_attrib("J kg-1", "cloud core specific internal energy", ""),
    )
    Variables["var_u_core"] = DiagnosticVariable(
        "var_u_core",
        diagnostic_var_attrib("m^2 s^-2", "cloud core variance of x-velocity", ""),
    )
    Variables["var_v_core"] = DiagnosticVariable(
        "var_v_core",
        diagnostic_var_attrib("m^2 s^-2", "cloud core variance of y-velocity", ""),
    )
    Variables["var_w_core"] = DiagnosticVariable(
        "var_w_core",
        diagnostic_var_attrib("m^2 s^-2", "cloud core variance of z-velocity", ""),
    )
    Variables["var_qt_core"] = DiagnosticVariable(
        "var_qt_core",
        diagnostic_var_attrib(
            "kg^2 kg^-2",
            "cloud core variance of total specific humidity",
            "",
        ),
    )
    Variables["var_thl_core"] = DiagnosticVariable(
        "var_thl_core",
        diagnostic_var_attrib(
            "K^2",
            "cloud core variance of liquid-ice potential temperature",
            "",
        ),
    )
    Variables["var_ei_core"] = DiagnosticVariable(
        "var_ei_core",
        diagnostic_var_attrib(
            "J^2 kg^-2",
            "cloud core variance of specific internal energy",
            "",
        ),
    )
    Variables["cov_w_rho_core"] = DiagnosticVariable(
        "cov_w_rho_core",
        diagnostic_var_attrib(
            "kg m^-2 s^-1",
            "cloud core vertical eddy flux of density",
            "",
        ),
    )
    Variables["cov_w_qt_core"] = DiagnosticVariable(
        "cov_w_qt_core",
        diagnostic_var_attrib(
            "kg kg^-1 m s^-1",
            "cloud core vertical eddy flux of specific humidity",
            "",
        ),
    )
    Variables["cov_w_thl_core"] = DiagnosticVariable(
        "cov_w_thl_core",
        diagnostic_var_attrib(
            "K m s^-1",
            "cloud core vertical eddy flux of liquid-ice potential temperature",
            "",
        ),
    )
    Variables["cov_w_ei_core"] = DiagnosticVariable(
        "cov_w_ei_core",
        diagnostic_var_attrib(
            "J kg^-1 m^-1 s^-1",
            "cloud core vertical eddy flux of specific internal energy",
            "",
        ),
    )
    Variables["cov_qt_thl_core"] = DiagnosticVariable(
        "cov_qt_thl_core",
        diagnostic_var_attrib(
            "kg kg^-1 K",
            "cloud core covariance of total specific humidity and liquid-ice potential temperature",
            "",
        ),
    )
    Variables["cov_qt_ei_core"] = DiagnosticVariable(
        "cov_qt_ei_core",
        diagnostic_var_attrib(
            "kg kg^-1 J kg^-1",
            "cloud core covariance of total specific humidity and specific internal energy",
            "",
        ),
    )
    Variables["E_k"] = DiagnosticVariable(
        "E_k",
        diagnostic_var_attrib(
            "",
            "volumetrically-averaged dimensionless kinetic energy",
            "",
        ),
    )
    Variables["dE"] = DiagnosticVariable(
        "dE",
        diagnostic_var_attrib(
            "",
            "volumetrically-averaged kinetic energy dissipation",
            "",
        ),
    )
    Variables["mass_loss"] =
        DiagnosticVariable("mass_loss", diagnostic_var_attrib("", "", ""))
    Variables["energy_loss"] =
        DiagnosticVariable("energy_loss", diagnostic_var_attrib("", "", ""))
=#
