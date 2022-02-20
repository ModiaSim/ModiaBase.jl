"""
Main module of ModiaBase.

* Developers: Hilding Elmqvist, Mogram AB, Martin Otter, DLR  
* First version: December 2020
* License: MIT (expat)

"""
module ModiaBase

const path    = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.9.0-dev"
const Date    = "2022-02-13"

#println("\nImporting ModiaBase Version $Version ($Date)")

using Unitful
using StaticArrays


# append! as needed in EquationAndStateInfo.jl and in ModiaLang/src/CodeGeneration.jl
appendVariable!(v1::Vector{FloatType}, s::FloatType)               where {FloatType}   = push!(v1,s)
appendVariable!(v1::Vector{FloatType}, v2::Vector{FloatType})      where {FloatType}   = append!(v1,v2)
appendVariable!(v1::Vector{FloatType}, v2::SVector{N,FloatType})   where {N,FloatType} = append!(v1,v2)
appendVariable!(v1::Vector{FloatType}, v2::NTuple{ N,FloatType})   where {N,FloatType} = append!(v1,v2)
@inline function appendVariable!(v1::Vector{FloatType}, v2::NTuple{N,SVector{M,FloatType}}) where {N,M,FloatType} 
    @inbounds for e in v2
        appendVariable!(v1,e)   # dispatch can be performed at compile-time, because typeof(e) = SVector{M,FloatType}
    end
end
@inline function appendVariable!(v1::Vector{FloatType}, v2::Tuple) where {FloatType}
    @inbounds for e in v2
        appendVariable!(v1,e)   # dispatch is performed at run-time, because typeof(e) is not known at compile-time.
    end
end


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
