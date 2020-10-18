##### Multiphysics tendency types

#####
##### First order fluxes
#####

"""
    Subsidence{FT} <: NonConservative

Subsidence
 - Equations: ρu
"""
struct Subsidence{FT} <: Source
    D::FT
end


#####
##### Second order fluxes
#####

struct CreateClouds <: Source end

