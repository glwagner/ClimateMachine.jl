using ..Atmos
using ..ConfigTypes
using ..DiagnosticsMachine

@diagnostics_group(
    "AtmosLESDefault",
    AtmosLESConfigType,
    GridDG,
    Nothing,
    # simple horizontal averages
    "u",
    "v",
    "w",
    "avg_rho",
    "rho",
    "temp",
    "pres",
    "thd",
    "et",
    "ei",
    "ht",
    "hi",
    "w_ht_sgs",
    # moisture related
    "qt",
    "ql",
    "qi",
    "qv",
    "thv",
    "thl",
    "w_qt_sgs",
    # variances and co-variances
    "u_u",
    "v_v",
    "w_w",
    "w_w_w",
    "ei_ei",
    "w_u",
    "w_v",
    "w_rho",
    "w_thd",
    "w_ei",
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
)

#= TODO
function atmos_les_default_clouds(
    ::MoistureModel,
    thermo,
    idx,
    qc_gt_0_z,
    qc_gt_0_full,
    z,
    cld_top,
    cld_base,
)
    return cld_top, cld_base
end
function atmos_les_default_clouds(
    moist::Union{EquilMoist, NonEquilMoist},
    thermo,
    idx,
    qc_gt_0_z,
    qc_gt_0_full,
    z,
    cld_top,
    cld_base,
)
    if thermo.moisture.has_condensate
        FT = eltype(qc_gt_0_z)
        qc_gt_0_z[idx] = one(FT)
        qc_gt_0_full[idx] = one(FT)

        cld_top = max(cld_top, z)
        cld_base = min(cld_base, z)
    end
    return cld_top, cld_base
end

function atmos_les_default_init(dgngrp::DiagnosticsGroup, currtime)
    atmos = Settings.dg.balance_law
    FT = eltype(Settings.Q)
    mpicomm = Settings.mpicomm
    mpirank = MPI.Comm_rank(mpicomm)

    atmos_collect_onetime(mpicomm, Settings.dg, Settings.Q)

    if mpirank == 0
        dims = OrderedDict("z" => (AtmosCollected.zvals, Dict()))

        # set up the variables we're going to be writing
        vars = OrderedDict()
        varnames = map(
            s -> startswith(s, "moisture.") ? s[10:end] : s,
            flattenednames(vars_atmos_les_default_simple(atmos, FT)),
        )
        ho_varnames = map(
            s -> startswith(s, "moisture.") ? s[10:end] : s,
            flattenednames(vars_atmos_les_default_ho(atmos, FT)),
        )
        append!(varnames, ho_varnames)
        for varname in varnames
            var = Variables[varname]
            vars[varname] = (("z",), FT, var.attrib)
        end
        vars["cld_frac"] = (("z",), FT, Variables["cld_frac"].attrib)
        vars["cld_top"] = ((), FT, Variables["cld_top"].attrib)
        vars["cld_base"] = ((), FT, Variables["cld_base"].attrib)
        vars["cld_cover"] = ((), FT, Variables["cld_cover"].attrib)
        vars["lwp"] = ((), FT, Variables["lwp"].attrib)

        # create the output file
        dprefix = @sprintf(
            "%s_%s_%s",
            dgngrp.out_prefix,
            dgngrp.name,
            Settings.starttime,
        )
        dfilename = joinpath(Settings.output_dir, dprefix)
        init_data(dgngrp.writer, dfilename, dims, vars)
    end

    return nothing
end

"""
    atmos_les_default_collect(dgngrp, currtime)

Collect the various 'AtmosLESDefault' diagnostic variables for the
current timestep and write them into a file.
"""
function atmos_les_default_collect(dgngrp::DiagnosticsGroup, currtime)
    mpicomm = Settings.mpicomm
    dg = Settings.dg
    Q = Settings.Q
    mpirank = MPI.Comm_rank(mpicomm)
    bl = dg.balance_law
    grid = dg.grid
    topology = grid.topology
    N = polynomialorder(grid)
    Nq = N + 1
    Nqk = dimensionality(grid) == 2 ? 1 : Nq
    npoints = Nq * Nq * Nqk
    nrealelem = length(topology.realelems)
    nvertelem = topology.stacksize
    nhorzelem = div(nrealelem, nvertelem)

    # get needed arrays onto the CPU
    if array_device(Q) isa CPU
        state_data = Q.realdata
        aux_data = dg.state_auxiliary.realdata
        vgeo = grid.vgeo
        ω = grid.ω
        gradflux_data = dg.state_gradient_flux.realdata
    else
        state_data = Array(Q.realdata)
        aux_data = Array(dg.state_auxiliary.realdata)
        vgeo = Array(grid.vgeo)
        ω = Array(grid.ω)
        gradflux_data = Array(dg.state_gradient_flux.realdata)
    end
    FT = eltype(state_data)

    zvals = AtmosCollected.zvals
    MH_z = AtmosCollected.MH_z

    # Visit each node of the state variables array and:
    # - generate and store the thermo variables,
    # - accumulate the simple horizontal sums, and
    # - determine the cloud fraction, top and base
    #
    thermo_array =
        [zeros(FT, num_thermo(bl, FT)) for _ in 1:npoints, _ in 1:nrealelem]
    simple_sums = [
        zeros(FT, num_atmos_les_default_simple_vars(bl, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    # for LWP
    ρq_liq_z = [zero(FT) for _ in 1:(Nqk * nvertelem)]
    # for cld*
    qc_gt_0_z = [zeros(FT, (Nq * Nq * nhorzelem)) for _ in 1:(Nqk * nvertelem)]
    qc_gt_0_full = zeros(FT, (Nq * Nq * nhorzelem))
    # In honor of PyCLES!
    cld_top = FT(-100000)
    cld_base = FT(100000)
    @visitQ nhorzelem nvertelem Nqk Nq begin
        evk = Nqk * (ev - 1) + k

        state = extract_state(dg, state_data, ijk, e, Prognostic())
        gradflux = extract_state(dg, gradflux_data, ijk, e, GradientFlux())
        aux = extract_state(dg, aux_data, ijk, e, Auxiliary())
        MH = vgeo[ijk, grid.MHid, e]

        thermo = thermo_vars(bl, thermo_array[ijk, e])
        compute_thermo!(bl, state, aux, thermo)

        simple = atmos_les_default_simple_vars(bl, simple_sums[evk])
        atmos_les_default_simple_sums!(
            bl,
            state,
            gradflux,
            aux,
            thermo,
            currtime,
            MH,
            simple,
        )

        idx = (Nq * Nq * (eh - 1)) + (Nq * (j - 1)) + i
        cld_top, cld_base = atmos_les_default_clouds(
            bl.moisture,
            thermo,
            idx,
            qc_gt_0_z[evk],
            qc_gt_0_full,
            zvals[evk],
            cld_top,
            cld_base,
        )

        # FIXME properly
        if isa(bl.moisture, EquilMoist) || isa(bl.moisture, NonEquilMoist)
            # for LWP
            ρq_liq_z[evk] += MH * thermo.moisture.q_liq * state.ρ * state.ρ
        end
    end

    # reduce horizontal sums and cloud data across ranks and compute averages
    simple_avgs = [
        zeros(FT, num_atmos_les_default_simple_vars(bl, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    cld_frac = zeros(FT, Nqk * nvertelem)
    for evk in 1:(Nqk * nvertelem)
        MPI.Allreduce!(simple_sums[evk], simple_avgs[evk], +, mpicomm)
        simple_avgs[evk] .= simple_avgs[evk] ./ MH_z[evk]

        # FIXME properly
        if isa(bl.moisture, EquilMoist) || isa(bl.moisture, NonEquilMoist)
            tot_qc_gt_0_z = MPI.Reduce(sum(qc_gt_0_z[evk]), +, 0, mpicomm)
            tot_horz_z = MPI.Reduce(length(qc_gt_0_z[evk]), +, 0, mpicomm)
            if mpirank == 0
                cld_frac[evk] = tot_qc_gt_0_z / tot_horz_z
            end

            # for LWP
            tot_ρq_liq_z = MPI.Reduce(ρq_liq_z[evk], +, 0, mpicomm)
            if mpirank == 0
                ρq_liq_z[evk] = tot_ρq_liq_z / MH_z[evk]
            end
        end
    end
    # FIXME properly
    if isa(bl.moisture, EquilMoist) || isa(bl.moisture, NonEquilMoist)
        cld_top = MPI.Reduce(cld_top, max, 0, mpicomm)
        if cld_top == FT(-100000)
            cld_top = NaN
        end
        cld_base = MPI.Reduce(cld_base, min, 0, mpicomm)
        if cld_base == FT(100000)
            cld_base = NaN
        end
        tot_qc_gt_0_full = MPI.Reduce(sum(qc_gt_0_full), +, 0, mpicomm)
        tot_horz_full = MPI.Reduce(length(qc_gt_0_full), +, 0, mpicomm)
        cld_cover = zero(FT)
        if mpirank == 0
            cld_cover = tot_qc_gt_0_full / tot_horz_full
        end
    end

    # complete density averaging
    simple_varnames = map(
        s -> startswith(s, "moisture.") ? s[10:end] : s,
        flattenednames(vars_atmos_les_default_simple(bl, FT)),
    )
    for evk in 1:(Nqk * nvertelem)
        simple_ha = atmos_les_default_simple_vars(bl, simple_avgs[evk])
        avg_rho = simple_ha.avg_rho
        for vari in 1:length(simple_varnames)
            if simple_varnames[vari] != "avg_rho"
                simple_avgs[evk][vari] /= avg_rho
            end
        end

        # for LWP
        # FIXME properly
        if isa(bl.moisture, EquilMoist) || isa(bl.moisture, NonEquilMoist)
            ρq_liq_z[evk] /= avg_rho
        end
    end

    # compute LWP
    lwp = NaN
    if mpirank == 0
        JcV = reshape(
            view(vgeo, :, grid.JcVid, topology.realelems),
            Nq^2,
            Nqk,
            nvertelem,
            nhorzelem,
        )
        Mvert = (ω .* JcV[1, :, :, 1])[:]
        lwp = FT(sum(ρq_liq_z .* Mvert))
    end

    # compute the variances and covariances
    ho_sums = [
        zeros(FT, num_atmos_les_default_ho_vars(bl, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    @visitQ nhorzelem nvertelem Nqk Nq begin
        evk = Nqk * (ev - 1) + k

        state = extract_state(dg, state_data, ijk, e, Prognostic())
        thermo = thermo_vars(bl, thermo_array[ijk, e])
        MH = vgeo[ijk, grid.MHid, e]

        simple_ha = atmos_les_default_simple_vars(bl, simple_avgs[evk])
        ho = atmos_les_default_ho_vars(bl, ho_sums[evk])
        atmos_les_default_ho_sums!(bl, state, thermo, MH, simple_ha, ho)
    end

    # reduce across ranks and compute averages
    ho_avgs = [
        zeros(FT, num_atmos_les_default_ho_vars(bl, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    for evk in 1:(Nqk * nvertelem)
        MPI.Reduce!(ho_sums[evk], ho_avgs[evk], +, 0, mpicomm)
        if mpirank == 0
            ho_avgs[evk] .= ho_avgs[evk] ./ MH_z[evk]
        end
    end

    # complete density averaging and prepare output
    if mpirank == 0
        varvals = OrderedDict()
        for (vari, varname) in enumerate(simple_varnames)
            davg = zeros(FT, Nqk * nvertelem)
            for evk in 1:(Nqk * nvertelem)
                davg[evk] = simple_avgs[evk][vari]
            end
            varvals[varname] = davg
        end

        ho_varnames = map(
            s -> startswith(s, "moisture.") ? s[10:end] : s,
            flattenednames(vars_atmos_les_default_ho(bl, FT)),
        )
        for (vari, varname) in enumerate(ho_varnames)
            davg = zeros(FT, Nqk * nvertelem)
            for evk in 1:(Nqk * nvertelem)
                simple_ha = atmos_les_default_simple_vars(bl, simple_avgs[evk])
                avg_rho = simple_ha.avg_rho
                davg[evk] = ho_avgs[evk][vari] / avg_rho
            end
            varvals[varname] = davg
        end

        if isa(bl.moisture, EquilMoist) || isa(bl.moisture, NonEquilMoist)
            varvals["cld_frac"] = cld_frac
            varvals["cld_top"] = cld_top
            varvals["cld_base"] = cld_base
            varvals["cld_cover"] = cld_cover
            varvals["lwp"] = lwp
        end

        # write output
        append_data(dgngrp.writer, varvals, currtime)
    end

    MPI.Barrier(mpicomm)
    return nothing
end # function collect
=#
