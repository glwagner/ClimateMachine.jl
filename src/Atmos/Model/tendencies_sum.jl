##### Sum wrappers

#####
##### First order fluxes
#####

function Σfluxes(fluxes::Tuple, m, state, aux, t, ts, direction)
    return sum(ntuple(Val(length(fluxes))) do i
        _flux_(fluxes[i], m, state, aux, t, ts, direction)
    end)
end

#####
##### Second order fluxes
#####

function Σfluxes(fluxes::Tuple, m, state, aux, t, ts, diffusive, hyperdiffusive)
    return sum(
        ntuple(Val(length(fluxes))) do i
            _flux_(fluxes[i], m, state, aux, t, ts, diffusive, hyperdiffusive)
        end,
    )
end

#####
##### Non-conservative sources
#####

function Σsources(sources::Tuple, m, state, aux, t, ts, direction, diffusive)
    return sum(
        ntuple(Val(length(sources))) do i
            _source_(sources[i], m, state, aux, t, ts, direction, diffusive)
        end,
    )
end
