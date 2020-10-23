using ..Atmos
using ..ConfigTypes
using ..DiagnosticsMachine

@diagnostics_group(
    "AtmosGCMDefault",
    AtmosGCMConfigType,
    GridInterpolated,
    nothing,
    # simple horizontal averages
    "u",
    "v",
    "w",
    "rho",
    "temp",
    "pres",
    "thd",
    "et",
    "ei",
    "ht",
    "hi",
    #"vort", TODO
    # moisture related
    "qt",
    "ql",
    "qv",
    "qi",
    "thv",
    "thl",
)
