module TestDifferentiate

using ModiaBase
using Test


removeBlock(ex) = ex
function removeBlock(ex::Expr)
    if ex.head in [:block]
        ex.args[1]
    else
        Expr(ex.head, [removeBlock(arg) for arg in ex.args]...)
    end
end

function showDifferentiate(ex, timeInvariants=[])
    Base.remove_linenums!(ex)
    print(removeBlock(ex), "   :   ")
    der = ModiaBase.derivative(ex, timeInvariants)
    Base.remove_linenums!(der)
    der = string(removeBlock(der))
    println(der)
    string(der)
end

function testDifferentiate()
    println("\nTest differentiate")

    @test showDifferentiate(10) == "0"
    @test showDifferentiate(10.0) == "0"
    @test showDifferentiate(:time) == "1"
    @test showDifferentiate(:x) == "der(x)"
    @test showDifferentiate(:(der(x))) == "der(der(x))"
    println()

    @test showDifferentiate(:(x.y.z)) == "der(x.y.z)"
    @test showDifferentiate(:(x = y)) == "der(x) = der(y)"
    println()

    @test showDifferentiate(:(+ y)) == "+(der(y))"
    @test showDifferentiate(:(x + y)) == "der(x) + der(y)"
    @test showDifferentiate(:(x + y + z)) == "der(x) + der(y) + der(z)"
    @test showDifferentiate(:(- y)) == "-(der(y))"
    @test showDifferentiate(:(x - y)) == "der(x) - der(y)"
    @test showDifferentiate(:(x - y - z)) == "(der(x) - der(y)) - der(z)"

    @test showDifferentiate(:(x * y)) == "der(x) * y + x * der(y)"
    @test showDifferentiate(:(x * y * z)) == "der(x) * y * z + x * der(y) * z + x * y * der(z)"
    @test showDifferentiate(:(x / y)) == "(one(x) / y) * der(x) + -((x / y) / y) * der(y)"
#    @test showDifferentiate(:(x / y / z)) == "(one(x / y) / z) * ((one(x) / y) * der(x) + -((x / y) / y) * der(y)) + -(((x / y) / z) / z) * der(z)"
    @test showDifferentiate(:(x ^ 2)) == "(2 * x ^ (2 - 1)) * der(x)"
    @test showDifferentiate(:(x ^ y)) ==  "(y * x ^ (y - 1)) * der(x) + if x isa Real && x <= 0\n            Base.oftype(float(x), NaN)\n        else\n            x ^ y * log(x)\n        end * der(y)"

#    @test showDifferentiate(:(x ^ y ^ z)) == "(y ^ z * x ^ (y ^ z - 1)) * der(x) + (x ^ (y ^ z) * log(x)) * ((z * y ^ (z - 1)) * der(y) + (y ^ z * log(y)) * der(z))"

    @test showDifferentiate(:(sin(x))) == "cos(x) * der(x)"

    @test showDifferentiate(:([x, y])) == "[der(x), der(y)]"

    @test showDifferentiate(:(if x; y else z end)) == "if x\n    der(y)\nelse\n    der(z)\nend"
    @test showDifferentiate(:(if x; y elseif z; v else w end)) == "if x\n    der(y)\nelseif z\n    der(v)\nelse\n    der(w)\nend"

    @test showDifferentiate(:(1+2+3)) == "0"
    @test showDifferentiate(:(1-2-3)) == "0"
    @test showDifferentiate(:(x-2-3)) == "der(x)"
    @test showDifferentiate(:(1-x-3)) == "-(der(x))"

    @test showDifferentiate(:(2*x + y + z)) == "2 * der(x) + der(y) + der(z)"

    # Previous test suit
    der = showDifferentiate(:(x + 5 + z = w))
    @test der == "der(x) + der(z) = der(w)"
       
#    der = showDifferentiate(ModiaBase.derivative(:(x + 5 + z = w)))
#    @test der == "der(der(x)) + der(der(z)) = der(der(w))"
    
    der = showDifferentiate(Expr(:(=), Expr(:call, :+, :x), :w))
    @test der == "+(der(x)) = der(w)"
    
    der = showDifferentiate(:(2 + 3 = w))
    @test der == "0 = der(w)"
    
    der = showDifferentiate(Expr(:(=), Expr(:call, :-, :x), :w))
    @test der == "-(der(x)) = der(w)"
    
    der = showDifferentiate(:(x - 5 - z = w))
    @test der == "der(x) - der(z) = der(w)"

    der = showDifferentiate(:(5 - x - z = w))
    @test der == "-(der(x)) - der(z) = der(w)"

    der = showDifferentiate(:(5x = w))
    @test der == "5 * der(x) = der(w)"
      
    der = showDifferentiate(:(x * 5 * z = w))
    @test der == "der(x) * 5 * z + x * 5 * der(z) = der(w)"
    
    der = showDifferentiate(:(4 * 5 * 6 = w))
    @test der == "0 = der(w)"
    
    der = showDifferentiate(:(y = x/5))  
    @test der == "der(y) = (one(x) / 5) * der(x)"
    
    der = showDifferentiate(:(y = 5/y))  
    @test der == "der(y) = -((5 / y) / y) * der(y)"
  
    der = showDifferentiate(:(y = [1, x]))
    @test der == "der(y) = [0, der(x)]"
    
    der = showDifferentiate(:(y = [2x 3x; 4x 5x]))
    @test der == "der(y) = [2 * der(x) 3 * der(x); 4 * der(x) 5 * der(x)]"
    
    der = showDifferentiate(:(y = [2*x 3x; 4x 5x]*[1, x]))
    @test der == "der(y) = [2 * der(x) 3 * der(x); 4 * der(x) 5 * der(x)] * [1, x] + [2x 3x; 4x 5x] * [0, der(x)]"
    #=
    der = showDifferentiate(:(y = transpose(B) + B´))
    @test der == "der(y) = transpose(der(B)) + der(B´)"
    =# 
    der = showDifferentiate(:(y = x[5, 6]))  
    @test der == "der(y) = (der(x))[5, 6]"
    
    der = showDifferentiate(:(y = x[5:7]))  
    @test der == "der(y) = (der(x))[5:7]"
    
  #=    der = showDifferentiate(:(y = [x for x in z]))  
    @test der == "der(y) = [x for x = der(z)]"
  
    der = showDifferentiate(:(y = [x[i] for i in 1:5]))  
    @test der == "der(y) = [x[i] for i = nothing]"
  =#      
    der = showDifferentiate(:(y = sin(x)))
    @test der == "der(y) = cos(x) * der(x)"
    
    der = showDifferentiate(:(y = cos(x)))
    @test der == "der(y) = -(sin(x)) * der(x)"
    
    #=
    der = showDifferentiate(:(y = tan(x)))
    @test der == "der(y) = (1 / cos(x) ^ 2) * der(x)"
    =#
    der = showDifferentiate(:(y = exp(x)))
    @test der == "der(y) = exp(x) * der(x)"
  
    der = showDifferentiate(:(z = x^y))
    #@test der == "der(z) = (y * x ^ (y - 1)) * der(x) + (x ^ y * log(x)) * der(y)"
    @test der == "der(z) = (y * x ^ (y - 1)) * der(x) + if x isa Real && x <= 0\n                Base.oftype(float(x), NaN)\n            else\n                x ^ y * log(x)\n            end * der(y)"
    
    der = showDifferentiate(:(y = log(x)))
    @test der == "der(y) = inv(x) * der(x)"
  
    der = showDifferentiate(:(y = asin(x)))  
    @test der == "der(y) = inv(sqrt(1 - x ^ 2)) * der(x)"
  
    der = showDifferentiate(:(y = acos(x)))  
    @test der == "der(y) = -(inv(sqrt(1 - x ^ 2))) * der(x)"
  
    der = showDifferentiate(:(y = atan(x)))  
    @test der == "der(y) = inv(1 + x ^ 2) * der(x)"

    #=
    der = showDifferentiate(:(y = f(x, 5, z)))
    @test der == "der(y) = f_der_1(x, 5, z) * der(x) + f_der_3(x, 5, z) * der(z)"
    
    der = showDifferentiate(:(y = f(x, 5, g(z))))
    @test der == "der(y) = f_der_1(x, 5, g(z)) * der(x) + f_der_3(x, 5, g(z)) * (g_der_1(z) * der(z))"
    =#
    der = showDifferentiate(:(y = true ? x : y))  
      @test der == "der(y) = if true\n        der(x)\n    else\n        der(y)\n    end"

    #=
    der = showDifferentiate(:(y = if b; x elseif false; y else z end))  
    @test der == "der(y) = if b\n        der(x)\n    else\n        if false\n            der(y)\n        else \n            der(z)\n        end\n    end"
    =#  
    der = showDifferentiate(:(y = time))  
    @test der == "der(y) = 1"
    
    der = showDifferentiate(:(y = a*x), [:a])  
    @test der == "der(y) = a * der(x)"
  end
  
  @testset "Symbolic" begin
  
    @testset "Differentiate" begin
      testDifferentiate()
    end
  
  end
  
end