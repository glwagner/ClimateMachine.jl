# # Eddy Diffusivity - Mass Flux module

using ClimateMachine.Atmos: nodal_update_auxiliary_state!

using ClimateMachine.DGMethods: LocalGeometry

import ClimateMachine.BalanceLaws:
    vars_state,
    init_state_prognostic!,
    flux_first_order!,
    flux_second_order!,
    compute_gradient_argument!,
    compute_gradient_flux!

import ClimateMachine.TurbulenceConvection:
    init_aux_turbconv!,
    turbconv_nodal_update_auxiliary_state!,
    turbconv_boundary_state!,
    turbconv_normal_boundary_flux_second_order!

vars_state(::Updraft, ::Prognostic, FT) = @vars(ρa::FT)

vars_state(m::NTuple{N, Updraft}, st::Prognostic, FT) where {N} =
    Tuple{ntuple(i -> vars_state(m[i], st, FT), N)...}

vars_state(m::EDMF, st::Prognostic, FT) =
    @vars(updraft::vars_state(m.updraft, st, FT))

function init_aux_turbconv!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    aux::Vars,
    geom::LocalGeometry,
) where {FT}
end;


function turbconv_nodal_update_auxiliary_state!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}
end;

function compute_gradient_argument!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}
end;

function compute_gradient_flux!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}

end;

function turbconv_source!(
    m::AtmosModel{FT},
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
) where {FT}
end;

function flux_first_order!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) where {FT}

end;

function flux_second_order!(
    turbconv::EDMF{FT},
    m::AtmosModel{FT},
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
) where {FT}
end;

function turbconv_boundary_state!(
    nf,
    bc::EDMFBCs,
    m::AtmosModel{FT},
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    state_int::Vars,
    aux_int::Vars,
) where {FT}
end;

function turbconv_normal_boundary_flux_second_order!(
    nf,
    bc::EDMFBCs,
    m::AtmosModel{FT},
    fluxᵀn::Vars,
    n⁻,
    state⁻::Vars,
    diff⁻::Vars,
    hyperdiff⁻::Vars,
    aux⁻::Vars,
    state⁺::Vars,
    diff⁺::Vars,
    hyperdiff⁺::Vars,
    aux⁺::Vars,
    bctype,
    t,
    _...,
) where {FT}
end;
