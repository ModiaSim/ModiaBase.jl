"""
Symbolic simplifications of +, -, *, /, ^

* Developer: Hilding Elmqvist, Mogram AB  
* First version: December 2020
* License: MIT (expat)

If possible, the operation is performed. Special care about adding zero and multiplying with one.
"""
module Simplify

export add, sub, mult, divide, power, simplify

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
    elseif x == -1 && ! isUnit(y)
        sub(0, y)  
    else
        :($x * $y)
    end
end

function divide(x, y)
    if typeof(x) in [Float64, Int64] && typeof(y) in [Float64, Int64]
        x / y
    elseif x == 0 || x == 0.0
        0
    elseif y === 1 || y == 1.0
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

simplify(ex) = ex

function simplify(ex::Expr)
    args = [simplify(arg) for arg in ex.args]
    if ex.head == :call && ex.args[1] in [:+, :-, :*] && all([typeof(a) in [Float64, Int64] for a in args[2:end]])
        if ex.args[1] in [:+]
            return sum(args[2:end])
        elseif ex.args[1] in [:-]
            return args[2] - sum(args[3:end])
        else
            return prod(args[2:end])
        end
    elseif ex.head == :call && ex.args[1] == :* && any([a == 0 for a in args[2:end]])
        # x*0*y = 0
        return 0
    elseif ex.head == :call && ex.args[1] == :^ && length(args) == 3 && args[3] == 1
        # x^1 = x
        return args[2]
    elseif ex.head == :call && ex.args[1] == :^ && length(args) == 3 && args[3] == 0
        # x^0 = 1
        return 1
    else
        return Expr(ex.head, args...)
    end
end

function testSimplify()
    @show simplify(:(1+3))
    @show simplify(:(3-1))
    @show simplify(:(4*3))
    @show simplify(:(4*3+x))
    @show simplify(:(1+3+x))
    @show simplify(:(1+x+3))
    @show simplify(:(x^0))
    @show simplify(:(x^1))
    @show simplify(:(x^(3-2)))
end

end