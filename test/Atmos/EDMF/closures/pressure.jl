#### Pressure model kernels
include(joinpath("..","helper_funcs", "diagnose_environment.jl"))

function perturbation_pressure(
    ss::AtmosModel{FT, N},
    m::PressureModel,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
    i::Int,
) where {FT, N}

    # Alias convention:
    gm = state
    en = state
    up = state.turbulence.updraft
    up_a = aux.turbulence.updraft
    up_d = diffusive.turbulence.updraft

    ρinv = 1 / gm.ρ
    en_area = 1 - sum([up[j].ρa for j in 1:N]) * ρinv
    w_env     = environment_w(state, aux, N)
    w_up = up[i].ρau[3] / up[i].ρa

    nh_press_buoy = -up[i].ρa * up_a[i].buoyancy * m.α_b
    nh_pressure_adv = up[i].ρa * m.α_a * w_up * up_d[i].∇u[3]
    nh_pressure_drag =
        -up[i].ρa * m.α_d * (w_up - w_env) * abs(w_up - w_env) / FT(500) #max(up_a[i].updraft_top, FT(500)) # this parameter should be exposed in the model

    dpdz = nh_press_buoy + nh_pressure_adv + nh_pressure_drag
    dpdz_tke_i = (w_up - w_env) * dpdz

    return dpdz, dpdz_tke_i
end;
