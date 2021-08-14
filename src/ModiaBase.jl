"""
Main module of ModiaBase.

* Developers: Hilding Elmqvist, Mogram AB, Martin Otter, DLR  
* First version: December 2020
* License: MIT (expat)

"""
module ModiaBase

const path    = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.7.2-dev"
const Date    = "2021-04-19"

#println("\nImporting ModiaBase Version $Version ($Date)")

using Unitful


import MonteCarloMeasurements
import StaticArrays

# Just temporary: should be removed
Base.append!(v::Vector{T}, s::T) where {T<:MonteCarloMeasurements.Particles}       = Base.push!(v,s)
Base.append!(v::Vector{T}, s::T) where {T<:MonteCarloMeasurements.StaticParticles} = Base.push!(v,s)

# append! as needed in ModiaBase
appendResidual!(v1::Vector{T}, s::T)          where {T} = push!(v1,s)
appendResidual!(v1::Vector{T}, v2::Vector{T}) where {T} = append!(v1,v2)
appendResidual!(v1::Vector{T}, v2::StaticArrays.SVector{N,T}) where {T,N} = append!(v1,v2)

include("LinearIntegerEquations.jl")

include("BLTandPantelidesUtilities.jl")
using .BLTandPantelidesUtilities

include("BLTandPantelides.jl")
using .BLTandPantelides

include("Differentiate.jl")
using .Differentiate

include("Tearing.jl")

include("Simplify.jl")
using .Simplify

include("Symbolic.jl")
using .Symbolic

include("EquationAndStateInfo.jl")
include("StateSelection.jl")

end
