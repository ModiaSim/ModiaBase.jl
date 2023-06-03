"""
Main module of ModiaBase.

* Developers: Hilding Elmqvist, Mogram AB, Martin Otter, DLR
* First version: December 2020
* License: MIT (expat)

"""
module ModiaBase

const path    = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.11.1"
const Date    = "2023-06-03"

#println("\nImporting ModiaBase Version $Version ($Date)")

using  Unitful
using  StaticArrays

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

include("NonlinearEquations.jl")
using .NonlinearEquations

end
