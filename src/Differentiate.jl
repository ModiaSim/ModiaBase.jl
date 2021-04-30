"""
Module for differentiation of expressions with regards to time.

* Developer: Hilding Elmqvist, Mogram AB  
* First version: December 2020
* License: MIT (expat)
"""
module Differentiate

using Base.Meta: isexpr
using DiffRules
using Unitful

export derivative

"""
    der = derivative(ex, timeInvariants=[])
Form the derivative of the expressions `ex` with regards to time. Time derivative of variables are denoted `der(v)`.

* `ex`: Expression to differentiate
* `timeInvariants`: Vector of time invariant variables, i.e. with time derivative equal to zero.
* `return `der`: Time derivative of ex 
"""
derivative(ex, timeInvariants=[]) = 0#u"1/s"
derivative(s::Symbol, timeInvariants=[]) = if s == :time; 1 elseif s in timeInvariants; 0 else :(der($s)) end # 0u"1/s"
function derivative(ex::Expr, timeInvariants=[])
    if isexpr(ex, :call) && ex.args[1] == :der
        :(der($ex))
    elseif isexpr(ex, :call) 
        func = ex.args[1]
        arguments = ex.args[2:end]
        derArguments = [derivative(a, timeInvariants) for a in arguments]
        if func in [:+, :-] 
            ders = [d for d in derArguments if d != 0]
            if length(ders) == 0 # der(1+2+3) == 0, der(1-2-3) == 0
                0
            elseif func == :- && length(derArguments) > 1 && length(ders) == 1 && derArguments[1] != 0 # der(x-2-3) == der(x)
                derArguments[1]
            else
                Expr(:call, func, ders...)
            end
        elseif func in [:*]
            # der(e1 * e2 * e3 + ...) = der(e1)*e2*e3 + e1*der(e2)*e3 + e1*e2*der(e3) + ...
            diff = Expr(:call, :+)
            for i in 1:length(arguments)
                terms = []
                for j in 1:length(arguments)
                    term = if i == j; derivative(arguments[j], timeInvariants) else arguments[j] end
                    if term == 0  
                        terms = []
                        break
                    else
                        push!(terms, term)
                    end
                end
                if length(terms) > 0
                    product = Expr(:call, :*, terms...)
                    push!(diff.args, product)   
                end             
            end   
            if length(diff.args) == 1 # Only + operator
                diff = 0
            elseif length(diff.args) == 2 # Only one + term, remove +
                diff = diff.args[2]
            else
                diff
            end
        elseif length(arguments) <= 2
            d = DiffRules.diffrule(:Base, ex.args[1], ex.args[2:end]...)
            if length(arguments) > 1 
                sum = []
                for i in 1:length(arguments)
                    if d[i] == :(one($(arguments[i]))) 
                        push!(sum, :($(derArguments[i])))
                    elseif d[i] == :(-one($(arguments[i])))
                        push!(sum, :(-$(derArguments[i])))
                    elseif derArguments[i] != 0
                        push!(sum, :($(d[i]) * $(derArguments[i])))
                    end
                end
                if length(sum) > 1
                    Expr(:call, :+, sum...)
                else
                    sum[1]
                end
            else
                :($d * $(derArguments[1]))
            end
        end
    elseif isexpr(ex, :.)
        if ex in timeInvariants; 0 else :(der($ex)) end
    elseif isexpr(ex, :if) || isexpr(ex, :elseif)
        Expr(ex.head, ex.args[1], [derivative(e, timeInvariants) for e in ex.args[2:end]]...) # Don't differentiate condition
    elseif isexpr(ex, :ref)
        Expr(ex.head, derivative(ex.args[1], timeInvariants), ex.args[2:end]...) # Don't differentiate indices
    else
        # For example: =, vect, hcat
        Expr(ex.head, [derivative(e, timeInvariants) for e in ex.args]...)
    end
end

end
