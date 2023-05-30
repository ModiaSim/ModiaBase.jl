module Runtests

using Test


@testset "Test ModiaBase" begin
    include("TestSymbolic.jl")
    include("TestBLTandPantelides.jl")
    include("TestDifferentiate.jl")
    include("TestTearing.jl")
    include("TestLinearIntegerEquations.jl")
    include("TestNonlinearEquations.jl")
end

end
