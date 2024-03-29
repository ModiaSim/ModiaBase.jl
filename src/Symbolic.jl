"""
Structural and symbolic processing.

* Developer: Hilding Elmqvist, Mogram AB  
* First version: December 2020
* License: MIT (expat)

Examples of use can be found in TestSymbolic.jl

"""
module Symbolic

export removeBlock, removeQuoteNode, makeDerVar, append, prepend, Incidence, findIncidence!, linearFactor, solveEquation, 
    isLinear, getCoefficients, substitute, removeUnits, resetEventCounters, getEventCounters, substituteForEvents

using Base.Meta: isexpr
using OrderedCollections
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
removeBlock(arr::Array{Any,1}) = [removeBlock(a) for a in arr]
removeBlock(arr::Array{Expr,1}) = [removeBlock(a) for a in arr]
removeBlock(q::QuoteNode) = q 
    # if typeof(q.value) == Symbol; q.value else q end # begin print("QuoteNode: "); dump(q); q; end # q.value
removeBlock(d::OrderedDict) = OrderedDict{Symbol, Any}([(k,removeBlock(v)) for (k,v) in d])

function removeBlock(ex::Expr)
    if isexpr(ex, :block) && length(ex.args) == 2 && typeof(ex.args[1]) == LineNumberNode 
        ex.args[2]
    elseif isexpr(ex, :block) && length(ex.args) == 3 && typeof(ex.args[1]) == LineNumberNode && typeof(ex.args[2]) == LineNumberNode 
        ex.args[3]
    else
        Expr(ex.head, [removeBlock(arg) for arg in ex.args]...)
    end
end

removeQuoteNode(ex) = ex
removeQuoteNode(arr::Array{Any,1}) = [removeQuoteNode(a) for a in arr]
removeQuoteNode(arr::Array{Expr,1}) = [removeQuoteNode(a) for a in arr]
removeQuoteNode(q::QuoteNode) = q.value
    # if typeof(q.value) == Symbol; q.value else q end # begin print("QuoteNode: "); dump(q); q; end # q.value
removeQuoteNode(d::OrderedDict) = OrderedDict{Symbol, Any}([(k,removeQuoteNode(v)) for (k,v) in d])

function removeQuoteNode(ex::Expr)
    Expr(ex.head, [removeQuoteNode(arg) for arg in ex.args]...)
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

function prepend(ex::Expr, prefix)
    if typeof(prefix) == Array{Any,1} && length(prefix) == 0
        ex
    elseif ex.head == :. && ex.args[1] == :up
        ex.args[2].value
    elseif ex.head in [:call, :kw]
        if false #ex.head == :call && ex.args[1] == :der
            e = Symbol(ex)
            :($prefix.$e)
        else
            Expr(ex.head, ex.args[1], [prepend(arg, prefix) for arg in ex.args[2:end]]...)
        end
    elseif ex.head == :macrocall
        if length(ex.args) >= 1 && ex.args[1] == Symbol("@u_str")
            # Do not expand  units, such as u"s", because otherwise logCode=true results in wrong code, because show(u"N") is N and not u"N".
            ex
        else
            eval(ex)
        end    
    else
        Expr(ex.head, [prepend(arg, prefix) for arg in ex.args]...)
    end
end


nCrossingFunctions = 0
nAfter = 0
nClocks = 0
nSamples = 0
previousVars = []
preVars = []
holdVars = []

function resetEventCounters()
    global nCrossingFunctions
    global nAfter
    global nClocks
    global nSamples 
	global previousVars
    global preVars
    global holdVars
    nCrossingFunctions = 0
    nAfter = 0
    nClocks = 0
    nSamples = 0
	previousVars = []
    preVars = []
    holdVars = []
end

function getEventCounters()
    global nCrossingFunctions
    global nAfter
    global nClocks
    global nSamples
	global previousVars
    global preVars
    global holdVars
    return (nCrossingFunctions, nAfter, nClocks, nSamples, previousVars, preVars, holdVars)
end

substituteForEvents(ex) = ex

function substituteForEvents(ex::Expr)
    global nCrossingFunctions
    global nAfter
    global nClocks
    global nSamples 
	global previousVar
    global preVars
    global holdVars
    if ex.head in [:call, :kw]
        if ex.head == :call && ex.args[1] == :positive
            nCrossingFunctions += 1
            :(positive(instantiatedModel, $nCrossingFunctions, ustrip($(substituteForEvents(ex.args[2]))), $(string(substituteForEvents(ex.args[2]))), _leq_mode))
        elseif ex.head == :call && ex.args[1] == :Clock
            @assert 2<=length(ex.args)<=3 "The Clock function takes one or two arguments: $ex"
             nClocks += 1
             if length(ex.args) == 2
                :(Clock(ustrip($(substituteForEvents(ex.args[2]))), instantiatedModel, $nClocks))
             else
                :(Clock(ustrip($(substituteForEvents(ex.args[2]))), ustrip($(substituteForEvents(ex.args[3]))), instantiatedModel, $nClocks))
             end
        elseif ex.head == :call && ex.args[1] == :sample
            nSamples += 1
            :(sample($(substituteForEvents(ex.args[2])), $(substituteForEvents(ex.args[3])), instantiatedModel, $nSamples))
        elseif ex.head == :call && ex.args[1] == :pre
            if length(ex.args) == 2
                push!(preVars, ex.args[2])
                nPre = length(preVars)
                :(pre(instantiatedModel, $nPre))
            else
                error("The pre function takes one arguments: $ex")
            end
        elseif ex.head == :call && ex.args[1] == :previous
            if length(ex.args) == 3
                push!(previousVars, ex.args[2])
                nPrevious = length(previousVars)
                :(previous($(substituteForEvents(ex.args[3])), instantiatedModel, $nPrevious))
            else
                error("The previous function presently takes two arguments: $ex")
            end
        elseif ex.head == :call && ex.args[1] == :hold
            push!(holdVars, ex.args[2])
            nHold = length(holdVars)
            if length(ex.args) == 3
                :(hold($(substituteForEvents(ex.args[2])), $(substituteForEvents(ex.args[3])), instantiatedModel, $nHold))
            else
#                error("The hold function takes two or three arguments, hold(v, clock) or hold(expr, start, clock) : $ex")
                error("The hold function takes two arguments, hold(v, clock): $ex")
            end
        elseif ex.head == :call && ex.args[1] in [:initial, :terminal]
            if length(ex.args) == 1
                :($(ex.args[1])(instantiatedModel))
            else
                error("The $(ex.args[1]) function don't take any arguments: $ex")
            end
        elseif ex.head == :call && ex.args[1] == :after
            # after(instantiatedModel, nr, t, tAsString, leq_mode) 
            nAfter += 1
            :(after(instantiatedModel, $nAfter, ustrip($(substituteForEvents(ex.args[2]))), $(string(substituteForEvents(ex.args[2]))), _leq_mode))
        else
            Expr(ex.head, ex.args[1], [substituteForEvents(arg) for arg in ex.args[2:end]]...)
        end
    else
        Expr(ex.head, [substituteForEvents(arg) for arg in ex.args]...)
    end
end

Incidence = Union{Symbol, Expr}

"""
    findIncidence!(ex, incidence::Array{Incidence,1})

Traverses an expression and finds incidences of Symbols and der(...)

* `ex`: Expression or array of expressions
* `incidence`: array of incidences. New incidences of `ex` are pushed. 
"""
findIncidence!(ex, incidence::Array{Incidence,1}, includeX::Bool=true) = nothing
findIncidence!(s::Symbol, incidence::Array{Incidence,1}, includeX::Bool=true) = begin if ! (s in [:(:), :end]); push!(incidence, s) end end
findIncidence!(arr::Array{Any,1}, incidence::Array{Incidence,1}, includeX::Bool=true) = [findIncidence!(a, incidence, includeX) for a in arr]
findIncidence!(arr::Array{Expr,1}, incidence::Array{Incidence,1}, includeX::Bool=true) = [findIncidence!(a, incidence, includeX) for a in arr]
function findIncidence!(ex::Expr, incidence::Array{Incidence,1}, includeX::Bool=true)
    if ex.head == :macrocall && ex.args[1] == Symbol("@u_str")
        nothing  
    elseif isexpr(ex, :function)
        nothing  
    elseif ex.head == :call
        if ex.args[1] == :der
            push!(incidence, ex) # der(x)
            if includeX
                push!(incidence, ex.args[2]) # x
            end
        elseif ex.args[1] in [:pre, :previous]
            [findIncidence!(e, incidence, includeX) for e in ex.args[3:end]] # skip operator/function name and first argument
		else
            [findIncidence!(e, incidence, includeX) for e in ex.args[2:end]] # skip operator/function name
        end
    elseif ex.head == :.
        push!(incidence, ex)
#		if ex.args[2].value != :all
#			push!(incidence, ex.args[1])
#		end
    elseif ex.head == :generator
        vars = [v.args[1] for v in ex.args[2:end]]
        incid = Incidence[]
        [findIncidence!(e, incid, includeX) for e in ex.args]
        unique!(incid)
        setdiff!(incid, vars)
        push!(incidence, incid...)
    else
        # For example: =, vect, hcat, block, ref
        [findIncidence!(e, incidence, includeX) for e in ex.args]
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
linearFactor(ex::Symbol, x::Incidence) =  if ex == x; (0, 1, true) else (ex, 0, true) end
function linearFactor(ex::Expr, x::Incidence)
    if ex.head == :block
        linearFactor(ex.args[1], x)    
    elseif ex.head == :macrocall && ex.args[1] == Symbol("@u_str")
        (ex, 0, true)    
    elseif isexpr(ex, :call) && ex.args[1] == :der
        if ex == x; (0, 1, true) else (ex, 0, true) end
    elseif isexpr(ex, :call) && ex.args[1] in [:positive, :previous]
        (ex, 0, true)
    elseif isexpr(ex, :call)
        func = ex.args[1]
        if func in [:_DUPLICATEEQUATION, :_DUPLICATIMPLICITDEPENDENCY, :implicitDependency]
            return (ex, 0, false)
        end
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
function isLinear(equ::Expr, x::Incidence)
    (rest, factor, linear) = linearFactor(equ, x)
    (linear, typeof(factor) != Expr)
end


"""
    (incidence, coefficients, rest, linear) = getCoefficients(ex)

If `ex` is `linear` with regards to all `incidence` (Symbols and der(...)), the `coefficients` and `rest` are returned 

* `ex`: Expression
* `return (incidence, coefficients, rest, linear)` 
"""

function getCoefficients(ex::Expr)
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
            break
        end
        (rest, factor, linearRest) = linearFactor(rest, x)
        if !linearRest
            linear = false
            break
        end
        push!(coefficients, factor)
    end
    incidence, coefficients, rest, linear
end

substitute(substitutions, ex) = begin nex = get(substitutions, ex, ex); if nex != ex; nex = substitute(substitutions, nex) else nex end; nex end 

substitute(substitutions, ex::Vector{Symbol}) = [substitute(substitutions, e) for e in ex]

substitute(substitutions, ex::Vector{Expr}) = [substitute(substitutions, e) for e in ex]

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

    println("TEST getCoefficients")
    n = 10000
    @show n
    e = Expr(:call, :f, fill(:x, n)...)
#    @show e
    @time incidence, coefficients, rest, linear = getCoefficients(e)
#    @show incidence coefficients rest linear

    e = Expr(:call, :_DUPLICATEEQUATION, fill(:x, n)...)
#    @show e
    @time incidence, coefficients, rest, linear = getCoefficients(e)
#    @show incidence coefficients rest linear
end

# testSubstitutions()

end