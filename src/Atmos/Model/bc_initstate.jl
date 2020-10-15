using ..Mesh.Grids: _x1, _x2, _x3
"""
    InitStateBC

Set the value at the boundary to match the `init_state_prognostic!` function. This is
mainly useful for cases where the problem has an explicit solution.

# TODO: This should be fixed later once BCs are figured out (likely want
# different things here?)
"""
struct InitStateBC end
function atmos_boundary_state!(
    ::Union{NumericalFluxFirstOrder, NumericalFluxGradient},
    bc::InitStateBC,
    m::AtmosModel,
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
    # Here we create a dummy_vgeo object so we can use LocalGeometry
    FT = eltype(aux⁺)
    dummy_vgeo = MArray{Tuple{1, _x3, 1}, FT}(undef)
    dummy_vgeo[1, _x1] = aux⁺.coord[1]
    dummy_vgeo[1, _x2] = aux⁺.coord[2]
    dummy_vgeo[1, _x3] = aux⁺.coord[3]
    init_state_prognostic!(
        m,
        state⁺,
        aux⁺,
        LocalGeometry{0, 0}(dummy_vgeo, 1, 1),
        t,
    )
end

function atmos_normal_boundary_flux_second_order!(
    nf,
    bc::InitStateBC,
    atmos,
    fluxᵀn,
    n⁻,
    state⁻,
    diff⁻,
    hyperdiff⁻,
    aux⁻,
    state⁺,
    diff⁺,
    hyperdiff⁺,
    aux⁺,
    bctype,
    t,
    args...,
)

    normal_boundary_flux_second_order!(
        nf,
        atmos,
        fluxᵀn,
        n⁻,
        state⁻,
        diff⁻,
        hyperdiff⁻,
        aux⁻,
        state⁺,
        diff⁺,
        hyperdiff⁺,
        aux⁺,
        bc,
        t,
        args...,
    )

end


function boundary_state!(
    ::NumericalFluxSecondOrder,
    m::AtmosModel,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    bc::InitStateBC,
    t,
    args...,
)
    # Here we create a dummy_vgeo object so we can use LocalGeometry
    FT = eltype(aux⁺)
    dummy_vgeo = MArray{Tuple{1, _x3, 1}, FT}(undef)
    dummy_vgeo[1, _x1] = aux⁺.coord[1]
    dummy_vgeo[1, _x2] = aux⁺.coord[2]
    dummy_vgeo[1, _x3] = aux⁺.coord[3]
    init_state_prognostic!(
        m,
        state⁺,
        aux⁺,
        LocalGeometry{0, 0}(dummy_vgeo, 1, 1),
        t,
    )
end
