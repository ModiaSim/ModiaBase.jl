module TestSymbolic

using ModiaBase
using ModiaBase.BLTandPantelidesUtilities
using ModiaBase.Symbolic

using Test

function showFindIncidence(ex)
    Base.remove_linenums!(ex)
    ex = removeBlock(ex)
    print(rpad(string(ex), 40))
    incidence = Incidence[]
    findIncidence!(ex, incidence)
    println(incidence)
    incidence
end

function testFindIncidence()
    println("testFindIncidence")
    @test showFindIncidence(10) == []
    @test showFindIncidence(:(x)) == [:x]
    @test showFindIncidence(:(x*(y + z*sin(w)))) == [:x, :y, :z, :w]
    @test showFindIncidence(:(x, x.y, x.y.z, x.y[z], f(x.y))) == [:x, :(x.y), :(x.y.z), :(x.y), :z, :(x.y)]
    println()
end

function showLinearFactor(ex, x)
    Base.remove_linenums!(ex)
    ex = removeBlock(ex)
    print(rpad(string(ex), 20))
    print(rpad(string(x), 10))
    factor = linearFactor(ex, x)
    Base.remove_linenums!.(factor)
    factor = string(removeBlock(factor))
    println(factor)
    string(factor)
end

function testLinearFactors()
    println("testLinearFactors")
    @test showLinearFactor(:(10), :x) == "(10, 0, true)"
    @test showLinearFactor(:(20.0), :x) == "(20.0, 0, true)"
    @test showLinearFactor(:(x), :x) == "(0, 1, true)"
    @test showLinearFactor(:(x.y.z), :(x.y.z)) == "(0, 1, true)"
    @test showLinearFactor(:(y), :x) == "(:y, 0, true)"
    @test showLinearFactor(:(x + y), :x) == "(:y, 1, true)"
    @test showLinearFactor(:(x + y + x + y), :x) == "(:(y + y), 2, true)"

    @test showLinearFactor(:(x - y), :x) == "(:(-y), 1, true)"
    @test showLinearFactor(:(y - x), :x) == "(:y, -1, true)"
    @test showLinearFactor(:(-x + y), :x) == "(:y, -1, true)"
    @test showLinearFactor(:(x - y - z), :x) == "(:(-y - z), 1, true)"
    @test showLinearFactor(:(x - y - x), :x) == "(:(-y), 0, true)"
    @test showLinearFactor(:(x - y + x), :x) == "(:(-y), 2, true)"

    @test showLinearFactor(:(x * y), :x) == "(0, :y, true)"
    @test showLinearFactor(:(x * x), :x) == "(0, 0, false)"
    @test showLinearFactor(:(2 * x), :x) == "(0, 2, true)"
    @test showLinearFactor(:(x * y * z), :x) == "(0, :(z * y), true)"

    @test showLinearFactor(:(x / y), :x) == "(0, :(1 / y), true)"
    @test showLinearFactor(:(y / x), :x) == "(:(y / 0), NaN, false)"
    @test showLinearFactor(:(x / y / z), :x) == "(0, :((1 / y) / z), true)"

    @test showLinearFactor(:(y \ x), :x) == "(0, :(1 / y), true)"

    @test showLinearFactor(:(sin(x)), :x) == "(:(sin(x)), 0, false)"
    @test showLinearFactor(:(sin(x) + x), :x) == "(:(sin(x)), 1, false)"
    @test showLinearFactor(:(sin(y)), :x) == "(:(sin(y)), 0, true)"

    @test showLinearFactor(:(if cond; x else y end), :x) == "(:(if cond\n      0\n  else\n      y\n  end), :(if cond\n      1\n  else\n      0\n  end), true)"
    @test showLinearFactor(:(if cond; x elseif cond2; y else z end), :x) ==  "(:(if cond\n      0\n  else\n      if cond2\n          y\n      else\n          z\n      end\n  end), :(if cond\n      1\n  else\n      0\n  end), true)"
    @test showLinearFactor(:(if cond; x+1 else 2x+2 end), :x) == "(:(if cond\n      1\n  else\n      2\n  end), :(if cond\n      1\n  else\n      2\n  end), true)"

    @test showLinearFactor(:(x[5]), :x) == "(:(x[5]), 0, false)"
    @test showLinearFactor(:(x = y), :x) == "(:(-y), 1, true)"
    @test showLinearFactor(:(x + 1 = y + 2), :x) == "(:(1 - (y + 2)), 1, true)"

    @test showLinearFactor(:(der(x) + x), :(der(x))) == "(:x, 1, true)"
    @test showLinearFactor(:(R*i = u), :i) == "(:(-u), :R, true)"
    println()
end

function showSolveEquation(ex, x)
    Base.remove_linenums!(ex)
    ex = removeBlock(ex)
    print(rpad(string(ex), 20))
    print(rpad(string(x), 10))
    factor = solveEquation(ex, x)
    Base.remove_linenums!.(factor)
    factor = string(removeBlock(factor))
    println(factor)
    string(factor)
end

function testSolveEquations()
    println("testSolveEquations")
    @test showSolveEquation(:(R*i = u), :i) == "(:(i = u / R), true)"
    @test showSolveEquation(:(R*i = u), :u) == "(:(u = R * i), true)"
    @test showSolveEquation(:(R*i = u), :R) == "(:(R = u / i), true)"
    println()

end

function showGetCoefficients(ex)
    Base.remove_linenums!(ex)
    ex = removeBlock(ex)
    print(rpad(string(ex), 20))
    factor = getCoefficients(ex)
    Base.remove_linenums!(factor)
    factor = string(removeBlock(factor))
    println(factor)
    string(factor)
end

function testGetCoefficients()
    println("testGetCoefficients")
    @test showGetCoefficients(:(v1 = 0)) == "(Union{Expr, Symbol}[:v1], Any[1], 0, true)"
    @test showGetCoefficients(:(v1 = 10)) == "(Union{Expr, Symbol}[:v1], Any[1], -10, true)"
    @test showGetCoefficients(:(v1 = v2)) == "(Union{Expr, Symbol}[:v1, :v2], Any[1, -1], 0, true)"
    @test showGetCoefficients(:(v1 = v2 + sin(v3))) == "(Union{Expr, Symbol}[:v1, :v2, :v3], Any[1, -1], :(-(sin(v3))), false)"
    @test showGetCoefficients(:(v1 = -v2)) == "(Union{Expr, Symbol}[:v1, :v2], Any[1, 1], 0, true)"
    @test showGetCoefficients(:(R*i = u)) == "(Union{Expr, Symbol}[:R, :i, :u], Any[:i], :(-u), false)"
    println()

end

@testset "Symbolic" begin

    @testset "TestFindIncidence" begin
        testFindIncidence()
    end

    @testset "TestLinearFactors" begin
        testLinearFactors()
    end

    @testset "TestSolveEquations" begin
        testSolveEquations()
    end

    @testset "TestGetCoefficients" begin
        testGetCoefficients()
    end
end


end