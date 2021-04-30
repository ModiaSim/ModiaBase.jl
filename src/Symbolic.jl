"""
Structural and symbolic processing.

* Developer: Hilding Elmqvist, Mogram AB  
* First version: December 2020
* License: MIT (expat)

Examples of use can be found in TestSymbolic.jl

"""
module Symbolic

export removeBlock, makeDerVar, append, prepend, Incidence, findIncidence!, linearFactor, solveEquation, 
    isLinear, getCoefficients, substitute, removeUnits, resetCounters

using Base.Meta: isexpr
#using OrderedCollections
using ModiaBase.Simplify
using Unitful
using Measurements
using MonteCarloMeasurements


"""
    e = removeBlock(ex)

Remove :block with LineNumberNode from expressions

* `ex`: Expression or array of expressions
* `return `e`: ex with :block removed 
"""
removeBlock(ex) = ex
#removeBlock(arr::Array{Any,1}) = [removeBlock(a) for a in arr]
removeBlock(arr::Array{Expr,1}) = [removeBlock(a) for a in arr]
function removeBlock(ex::Expr)
    if isexpr(ex, :block) && typeof(ex.args[1]) == LineNumberNode && length(ex.args) == 2
        ex.args[2]
    else
        Expr(ex.head, [removeBlock(arg) for arg in ex.args]...)
    end
end

"""
    e = makeDerVar(ex)

Recursively converts der(x) to Symbol(:(der(x))) in expression `ex`

* `ex`: Expression or array of expressions
* `return `e`: ex with der(x) converted 
"""
makeDerVar(ex, parameters) = if typeof(ex) in [Symbol, Expr] && ex in parameters; prepend(ex, :(_p)) else ex end

function makeDerVar(ex::Expr, parameters=[])
    if ex.head == :call && ex.args[1] == :der
        Symbol(ex)
	elseif isexpr(ex, :.) && ex in parameters
		prepend(ex, :(_p))
    elseif ex.head == :.
        Symbol(ex)
    else
        Expr(ex.head, [makeDerVar(arg, parameters) for arg in ex.args]...)
    end
end

removeUnits(ex) = if typeof(ex) <: Unitful.Quantity; @show ex; ustrip(ex) else ex end
removeUnits(ex::Expr) = if ex.head == :macrocall && ex.args[1] == Symbol("@u_str"); 1 else Expr(ex.head, [removeUnits(arg) for arg in ex.args]...) end

append(ex, suffix) = if ex == nothing; suffix else Expr(:., ex, QuoteNode(suffix)) end

prepend(ex, prefix) = ex
#prepend(ex, prefix::Array{Any,1}) = if length(prefix) == 0; ex else prepend(prepend(ex, prefix[end]), prefix[1:end-1]) end
#prepend(ex, prefix::Array{Symbol,1}) = if length(prefix) == 0; ex else prepend(prepend(ex, prefix[end]), prefix[1:end-1]) end
#prepend(ex::Symbol, prefix::Array{Symbol,1}) = if length(prefix) == 0; ex else prepend(prepend(ex, prefix[end]), prefix[1:end-1]) end
#prepend(ex::QuoteNode, prefix::Symbol) = Expr(:., prefix, QuoteNode(ex))
prepend(ex::Symbol, prefix) = if prefix == nothing; ex elseif ex in [:time, :instantiatedModel, :_leq_mode, :_x]; ex else Expr(:., prefix, QuoteNode(ex)) end
prepend(ex::Symbol, prefix::Nothing) = ex
prepend(arr::Array{Expr,1}, prefix) = [prepend(a, prefix) for a in arr]
#prepend(dict::OrderedCollections.OrderedDict{Symbol,Expr}, prefix) = OrderedDict([prepend(k, prefix) => prepend(k, prefix) for (k,v) in dict])

nCrossingFunctions = 0
nClocks = 0
nSamples = 0

function resetCounters()
    global nCrossingFunctions
    global nClocks
    global nSamples 
    nCrossingFunctions = 0
    nClocks = 0
    nSamples = 0
end

function prepend(ex::Expr, prefix)
    global nCrossingFunctions
    global nClocks
    global nSamples 
    if typeof(prefix) == Array{Any,1} && length(prefix) == 0
        ex
    elseif ex.head == :. && ex.args[1] == :up
        ex.args[2].value
    elseif ex.head in [:call, :kw]
        if false #ex.head == :call && ex.args[1] == :der
            e = Symbol(ex)
            :($prefix.$e)
        elseif ex.head == :call && ex.args[1] == :positive
            nCrossingFunctions += 1
            :(positive(instantiatedModel, $nCrossingFunctions, $(prepend(ex.args[2], prefix)), $(string(prepend(ex.args[2], prefix))), _leq_mode))
        elseif ex.head == :call && ex.args[1] == :Clock
            nClocks += 1
            :(Clock($(prepend(ex.args[2], prefix)), instantiatedModel, $nClocks))
        elseif ex.head == :call && ex.args[1] == :sample
            nSamples += 1
            :(sample($(prepend(ex.args[2], prefix)), $(prepend(ex.args[3], prefix)), instantiatedModel, $nSamples))
        else
            Expr(ex.head, ex.args[1], [prepend(arg, prefix) for arg in ex.args[2:end]]...)
        end
    elseif ex.head == :macrocall
        eval(ex)
    else
        Expr(ex.head, [prepend(arg, prefix) for arg in ex.args]...)
    end
end


Incidence = Union{Symbol, Expr}

"""
    findIncidence!(ex, incidence::Array{Incidence,1})

Traverses an expression and finds incidences of Symbols and der(...)

* `ex`: Expression or array of expressions
* `incidence`: array of incidences. New incidences of `ex` are pushed. 
"""
findIncidence!(ex, incidence::Array{Incidence,1}) = nothing
findIncidence!(s::Symbol, incidence::Array{Incidence,1}) = begin if s != :(:); push!(incidence, s) end end
findIncidence!(arr::Array{Any,1}, incidence::Array{Incidence,1}) = [findIncidence!(a, incidence) for a in arr]
findIncidence!(arr::Array{Expr,1}, incidence::Array{Incidence,1}) = [findIncidence!(a, incidence) for a in arr]
function findIncidence!(ex::Expr, incidence::Array{Incidence,1})
    if ex.head == :macrocall && ex.args[1] == Symbol("@u_str")
        nothing    
    elseif ex.head == :call
        if ex.args[1] == :der
            push!(incidence, ex) # der(x)
            push!(incidence, ex.args[2]) # x
        else
            [findIncidence!(e, incidence) for e in ex.args[2:end]] # skip operator/function name
        end
    elseif ex.head == :.
        push!(incidence, ex)
#		if ex.args[2].value != :all
#			push!(incidence, ex.args[1])
#		end
    elseif ex.head == :generator
        vars = [v.args[1] for v in ex.args[2:end]]
        incid = Incidence[]
        [findIncidence!(e, incid) for e in ex.args]
        unique!(incid)
        setdiff!(incid, vars)
        push!(incidence, incid...)
    else
        # For example: =, vect, hcat, block, ref
        [findIncidence!(e, incidence) for e in ex.args]
    end
    nothing
end


"""
    (rest, factor, linear) = linearFactor(ex, x)

Finds the linear `factor` and `rest` if `ex` is `linear` with regards to `x` (ex == factor*x + rest)

* `ex`: Expression 
* `return (rest, factor, linear)`:  
"""
linearFactor(ex, x) = (ex, 0, true)
linearFactor(ex::Symbol, x) =  if ex == x; (0, 1, true) else (ex, 0, true) end
function linearFactor(ex::Expr, x)
    if ex.head == :block
        linearFactor(ex.args[1], x)    
    elseif ex.head == :macrocall && ex.args[1] == Symbol("@u_str")
        (ex, 0, true)    
    elseif isexpr(ex, :call) && ex.args[1] == :der
        if ex == x; (0, 1, true) else (ex, 0, true) end
    elseif isexpr(ex, :call)
        func = ex.args[1]
        arguments = ex.args[2:end]
        factored = [linearFactor(a, x) for a in arguments]
        rests = [f[1] for f in factored]
        factors = [f[2] for f in factored]
        linears = [f[3] for f in factored]
        if func == :+
            rest = foldl(add, rests)
            factor = foldl(add, factors)
            (rest, factor, all(linears))
        elseif func == :-
            if length(rests) == 1
                rest = sub(0, rests[1])
            else
                rest = foldl(sub, rests)
            end
            if length(factors) == 1
                factor = sub(0, factors[1])
            else
                factor = foldl(sub, factors)
            end
            (rest, factor, all(linears))
        elseif func == :*
            if length(arguments) > 2
                linearFactor(foldl(mult, arguments), x)
            else
                # (r1 + f1*x)*(r2 + f2*x) = (r1*r2 + r1*f2*x + f1*r2*x + ...) 
                rest = foldl(mult, rests)
                factor = 0
                if factors[1] == 0
                    factor = mult(rests[1], factors[2])
                    (rest, factor, all(linears))
                elseif length(factors) == 2 && factors[2] == 0
                    factor = mult(rests[2], factors[1])
                    (rest, factor, all(linears))
                else
                    (0, 0, false)
                end
            end
        elseif func == :/
            # (r1 + f1*x)/r2 = (r1/r2 + f1/r2*x) 
            rest = foldl(divide, rests)
            @assert length(factors) == 2 "Non-binary division is not handled."
            factor = divide(factors[1], rests[2])
            (rest, factor, all(linears) && all(factors[2:end] .== 0))
        elseif func == :\
            # r1 \ (r2 + f2*x) = (r1\r2 + r1\f2*x) 
            rest = divide(rests[2], rests[1])
            @assert length(factors) == 2 "Non-binary \\ is not handled."
            factor = divide(factors[2], rests[1])
            (rest, factor, all(linears) && factors[1] == 0)
        else
            (ex, 0, all(linears) && all(factors .== 0))
        end
    elseif ex.head == :.
        if ex == x; (0, 1, true) else (ex, 0, true) end
    elseif isexpr(ex, :if) || isexpr(ex, :elseif)
        cond = linearFactor(ex.args[1], x)
        then = linearFactor(ex.args[2], x)
        els = linearFactor(ex.args[3], x)
        if then[1] == els[1]
            rest = then[1]
        else
            rest = :(if $(cond[1]); $(then[1]) else $(els[1]) end)
        end
        if then[2] == els[2]
            factor = then[2]
        else
            factor = :(if $(cond[1]); $(then[2]) else $(els[2]) end)
        end
        (rest, factor, cond[3] && then[3] && els[3])
    elseif isexpr(ex, :(=))
        LHS = linearFactor(ex.args[1], x)
        RHS = linearFactor(ex.args[2], x)
        rest = sub(LHS[1], RHS[1])
        factor = sub(LHS[2], RHS[2])
        (rest, factor, LHS[3] && RHS[3])
    elseif isexpr(ex, :vect) || isexpr(ex, :vcat) || isexpr(ex, :hcat) || isexpr(ex, :row)
        arguments = ex.args[2:end]
        factored = [linearFactor(a, x) for a in arguments]
        linears = [f[3] for f in factored]
        (ex, 0, all(linears))
    else
#        @warn "Unknown expression type" ex
#        dump(ex)
        incidence = Incidence[]
        findIncidence!(ex, incidence)
        if x in incidence
            (ex, 0, false)
        else
            (ex, 0, true)
        end
    end
end


"""
    (solved_equation, solved) = solveEquation(equ::Expr, x)

Solves `equ` for `x` if `ex` is linear with regards to `x` (equ == factor*x + rest = 0)

* `equ`: Equation (expression = expression)
* `return (solved_equation, solved)`: If the equation is `solved`, solved_equation constains the corresponding Expr 
"""
function solveEquation(equ::Expr, x)
    (rest, factor, linear) = linearFactor(equ, x)
    (:($x = $(divide(sub(0, rest), factor))), linear)
end


"""
    (linear, constant) = isLinear(ex::Expr, x)

Tests if `ex` is `linear` with regards to `x` (ex == factor*x + rest) and checks if the factor is `constant`

* `ex`: Expression
* `return (linear, constant)` 
"""
function isLinear(equ::Expr, x)
    (rest, factor, linear) = linearFactor(equ, x)
    (linear, typeof(factor) != Expr)
end


"""
    (incidence, coefficients, rest, linear) = getCoefficients(ex)

If `ex` is `linear` with regards to all `incidence` (Symbols and der(...)), the `coefficients` and `rest` are returned 

* `ex`: Expression
* `return (incidence, coefficients, rest, linear)` 
"""

function getCoefficients(ex)
    incidence = Incidence[]
    findIncidence!(ex, incidence)
    coefficients = []
    rest = ex
    linear = true
    for x in incidence
        crossIncidence = Incidence[]
        findIncidence!(coefficients, crossIncidence)
        if x in crossIncidence
            linear = false
        end
        (rest, factor, linearRest) = linearFactor(rest, x)
        if !linearRest
            linear = false
        end
        push!(coefficients, factor)
    end
    incidence, coefficients, rest, linear
end

substitute(substitutions, ex) = begin nex = get(substitutions, ex, ex); if nex != ex; nex = substitute(substitutions, nex) else nex end; nex end 

substitute(substitutions, ex::MonteCarloMeasurements.StaticParticles) = ex

substitute(substitutions, ex::Measurements.Measurement{Float64}) = ex

function substitute(substitutions, ex::Expr)
    if isexpr(ex, :quote)
        ex
    elseif ex.head == :.
        nex = get(substitutions, ex, ex)
        if nex != ex
            nex = substitute(substitutions, nex)
        end
        nex
    elseif ex.head == :call && ex.args[1] == :+ && length(ex.args) ==3
        add(substitute(substitutions, ex.args[2]), substitute(substitutions,ex.args[3]))
    elseif ex.head == :call && ex.args[1] == :- && length(ex.args) ==3
        sub(substitute(substitutions, ex.args[2]), substitute(substitutions,ex.args[3]))
    elseif ex.head == :call && ex.args[1] == :* && length(ex.args) ==3
        mult(substitute(substitutions, ex.args[2]), substitute(substitutions,ex.args[3]))
    else
        Expr(ex.head, [substitute(substitutions, arg) for arg in ex.args]...)
    end
end

function testSubstitutions()
    ex = substitute(Dict(:a=>:b), :(a+b=0))
    @show ex
    ex = substitute(Dict(:a=>:b), :(a+b+c=0))
    @show ex
    ex = substitute(Dict(:b=>0), :(a+b=0))
    @show ex
    ex = substitute(Dict(:b=>0.0), :(a+b=0))
    @show ex
    ex = substitute(Dict(:b=>0), :(a+b+c=0))
    @show ex
    ex = substitute(Dict(:b=>0.0), :(a+b+c=0))
    @show ex
    ex = substitute(Dict(:b=>0.0), :(a+b*(c+d)=0))
    @show ex
end

# testSubstitutions()

end