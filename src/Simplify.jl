"""
Symbolic simplifications of +, -, *, /, ^

* Developer: Hilding Elmqvist, Mogram AB  
* First version: December 2020
* License: MIT (expat)

If possible, the operation is performed. Special care about adding zero and multiplying with one.
"""
module Simplify

export add, sub, mult, divide, power

using Base.Meta: isexpr
using Unitful

isUnit(x) = typeof(x) <: Unitful.FreeUnits || isexpr(x, :macrocall)

function add(x, y)
    if typeof(x) in [Float64, Int64] && typeof(y) in [Float64, Int64]
        x + y
    elseif x == 0
        y
    elseif y == 0 
        x
    else
        :($x + $y)
    end
end

function sub(x, y)
    if typeof(x) in [Float64, Int64] && typeof(y) in [Float64, Int64]
        x - y
    elseif x == 0 && isexpr(y, :call) && length(y.args) == 2 && y.args[1] == :-
        y.args[2]
    elseif x == 0
        :(- $y)
    elseif y == 0 
        x
    else
        :($x - $y)
    end
end

function mult(x, y)
    if typeof(x) in [Float64, Int64] && typeof(y) in [Float64, Int64]
        x * y
    elseif (x == 0 && ! isUnit(y)) || (y == 0 && ! isUnit(x))  # Does not remove unit even if 0
        0
    elseif x == 1 && ! isUnit(y)
        y
    elseif y == 1 && ! isUnit(x)
        x  
    else
        :($x * $y)
    end
end

function divide(x, y)
    if typeof(x) in [Float64, Int64] && typeof(y) in [Float64, Int64]
        x / y
    elseif x === 0
        0
    elseif y === 1 
        x  
    elseif y === -1
        sub(0, x)  
    else
        :($x / $y)
    end
end

function power(x, y)
    if typeof(x) in [Float64, Int64] && typeof(y) in [Float64, Int64]
        x ^ y
    elseif y == 1 
        x
    else
        :($x ^ $y)
    end
end

end