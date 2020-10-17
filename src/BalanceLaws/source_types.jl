#### Source types

export AbstractSource, Flux1ˢᵗOrder, Flux2ⁿᵈOrder, NonConservative

"""
    AbstractSource

Subtypes of this describe tendencies, which
are partitioned into hierarchies for dispatching.
"""
abstract type AbstractSource end

abstract type Flux1ˢᵗOrder <: AbstractSource end
abstract type Flux2ⁿᵈOrder <: AbstractSource end
abstract type NonConservative <: AbstractSource end
