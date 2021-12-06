"""
Main module of ModiaBase.

* Developers: Hilding Elmqvist, Mogram AB, Martin Otter, DLR  
* First version: December 2020
* License: MIT (expat)

"""
module ModiaBase

const path    = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.8.0-dev"
const Date    = "2021-12-04"

#println("\nImporting ModiaBase Version $Version ($Date)")

using Unitful
import StaticArrays


# append! as needed in EquationAndStateInfo.jl
appendResidual!(v1::AbstractVector, s::Number)          = push!(v1,s)
appendResidual!(v1::AbstractVector, v2::AbstractVector) = append!(v1,v2)

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
