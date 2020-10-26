#### Source types

export AbstractTendency, Flux1ˢᵗOrder, Flux2ⁿᵈOrder, NonConservative
export AbstractTendencyDefinition,
    Flux1ˢᵗOrderTendency, Flux2ⁿᵈOrderTendency, NonConservativeTendency

"""
    AbstractTendency

Subtypes of this are used for specifying
a tuple of tendencies to be accumulated.
"""
abstract type AbstractTendency end
struct Flux1ˢᵗOrder <: AbstractTendency end
struct Flux2ⁿᵈOrder <: AbstractTendency end
struct NonConservative <: AbstractTendency end

"""
    AbstractTendencyDefinition

Subtypes of this are used for subtyping
each tendency definition.
"""
abstract type AbstractTendencyDefinition end
abstract type Flux1ˢᵗOrderTendency <: AbstractTendencyDefinition end
abstract type Flux2ⁿᵈOrderTendency <: AbstractTendencyDefinition end
abstract type NonConservativeTendency <: AbstractTendencyDefinition end
