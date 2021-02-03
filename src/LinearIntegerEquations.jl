# License for this file: MIT (expat)
# Copyright 2020-2021, DLR Institute of System Dynamics and Control
# Author: Martin Otter, DLR-SR

#=
Algorithm to exactly process linear Integer equations and hereby removing underdeterminism, overdeterminism,
and state constraints that are not structurally visible, as well as eliminating simple equations. 

The following is exported:

- [`simplifyLinearIntegerEquations`](@ref) 
- [`printSimplifiedLinearIntegerEquations`](@ref)

and the following utility functions that can be applied on the return
argument `vProperty` of [`simplifyLinearIntegerEquations`](@ref):

- `isNotEliminated(vProperty, v)` - if variable v is not eliminated.
- `isEliminated(vProperty, v)` - if variable v is eliminated.
- `isZero(vProperty,v)` - if eliminated variable v is zero.
- `isAlias(vProperty,v)` - if eliminated variable v is an alias variable v = v_alias
- `isNegAlias(vProperty,v)` - if eliminated variable v is a negative alias variable v = -v_alias
- `alias(vProperty,v)` - alias variable v_alias of eliminated variable v (v = v_alias).
- `negAlias(vProperty,v)` - negated alias variable v_alias of eliminated variable v (v = -v_alias).


# Main developer

[Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/), 
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
=#

import DataStructures

export simplifyLinearIntegerEquations!
export printSimplifiedLinearIntegerEquations

export isNotEliminated, isEliminated, isZero, isAlias, isNegAlias, alias, negAlias


# Inquire properties of eliminated variables (vEliminated is returned from simplifyLinearIntegerEquations(..)
const IS_PRESENT = typemin(Int)   # Variable is not eliminated


isNotEliminated(vProperty::Vector{Int}, v::Int) =  vProperty[v] == IS_PRESENT
isEliminated(   vProperty::Vector{Int}, v::Int) =  vProperty[v] != IS_PRESENT
isZero(         vProperty::Vector{Int}, v::Int) =  vProperty[v] == 0
isAlias(        vProperty::Vector{Int}, v::Int) =  vProperty[v] > 0
isNegAlias(     vProperty::Vector{Int}, v::Int) =  vProperty[v] < 0 && vProperty[v] != IS_PRESENT
alias(          vProperty::Vector{Int}, v::Int) =  vProperty[v]
negAlias(       vProperty::Vector{Int}, v::Int) = -vProperty[v]
    


"""
    getFirstActiveIndex(Gint, vActive, i::Int)
    
Return the first index j of equation i, such that vActive[ Gint[i][j] ] == true.
If there is no such index, return zero.
"""
function getFirstActiveIndex(Gint, vActive::Vector{Bool}, i::Int)::Int
    found = false
    j = 0
    for (k,vk) in enumerate(Gint[i])
        if vActive[vk]
            j = k
            break
        end
    end
    return j
end



""" 
    swapEquations!(Gint, eInt, GcInt, GcInt2, i, j)
    
Swap equations i and j
"""
function swapEquations!(Gint, eInt::Vector{Int}, GcInt, GcInt2, i::Int, j::Int)::Nothing
    Gint_i   = Gint[i]
    GcInt_i  = GcInt[i]    
    GcInt2_i = GcInt2[i]
    eInt_i   = eInt[i]
    
    Gint[i]   = Gint[j]
    GcInt[i]  = GcInt[j]    
    GcInt2[i] = GcInt2[j]
    eInt[i]   = eInt[j]
    
    Gint[j]   = Gint_i
    GcInt[j]  = GcInt_i    
    GcInt2[j] = GcInt2_i
    eInt[j]   = eInt_i
    return nothing
end


"""
    cj = find_v(eq_i, eqc_i, vj)
    
Return coefficient of variable vj in eq_i or zero, if vj not in eq_i.
"""
function find_v(eq_i, eqc_i, vj)::Int
    for (j, vi) in enumerate(eq_i)
        if vi == vj
            return eqc_i[j]
        end
    end
    return 0
end
        
mutable struct Buffer
    Gint_i::Vector{Int}
    GcInt2_i::Vector{Int}
    
    Buffer() = new(Int[], Int[])
end

""" 
    eliminateVariable!(Gint, GcInt2, k, k_perm, i, pivot, vj, oldPivot)

Eliminate variable vj from equation i and store the result as equation i.

k_perm is the sorting order of equation k (return argument of sortperm(..)).

If equation i does not contain variable vj, the equation is not changed.

If equation i contains variable vj, then multiply equation i with pivot and
equation k with the coefficient of variable vj in equation i.
Subtract the latter from the former equation and divide by oldPivot
(it is guaranteed that the division has no reminder).

If Gc[i,j] would be the coefficient of variable j in equation i, then the
following operation is carried out:

pivot_i = Gc[i,vj]
for j = 1:nv  # nv: number of variables
    Gc[i,j] = div(pivot * Gc[i,j] - pivot_i * Gc[k,j], oldPivot)  # guaranteed no remainder
end
        
If one of the coeffients of equation i, so GcInt2[i], becomes zero, then this coefficient and the corresponding
variable is removed from the equation. It is therefore guaranteed that after the operation,
no coefficient of the new equation i is zero, so no element of GcInt2[i][:] is zero
and no element of Gint[i][:] is vj.
"""
function eliminateVariable!(Gint, GcInt2, k::Int, i::Int, pivot::Int, vj::Int, oldPivot::Int, buffer::Buffer)::Nothing
    # Return if equation i does not contain variable vj   
    vj_present=false
    c_vj = 0   # Coefficient of vj
    for (index, value) in enumerate(Gint[i])
        if value == vj
            vj_present = true
            c_vj = GcInt2[i][index]
            deleteat!(Gint[i] ,index)
            deleteat!(GcInt2[i],index) 
            break
        end
    end
    if !vj_present
        return nothing
    end

    #println("... Eliminate variable $vj from equation $i = ", Gint[i])
    #println("    Gint[$k] = ", Gint[k], ", GcInt2[$k] = ", GcInt2[k], ", pivot = $pivot, c_vj = $c_vj")
    
    # Union of relevant variables
    eq_i   = Gint[i]
    eq_k   = Gint[k]
    eqc_i  = GcInt2[i]    
    eqc_k  = GcInt2[k]
    v_all  = union(eq_i, eq_k)
    #println("v_all = $v_all")

    # Eliminate variable vj from equation i
    empty!(buffer.Gint_i)
    empty!(buffer.GcInt2_i)    
    for v in v_all
        ci = find_v(eq_i, eqc_i, v)
        ck = find_v(eq_k, eqc_k, v)
        ci_new = div(pivot*ci - c_vj*ck, oldPivot) 
        if ci_new != 0
            push!(buffer.Gint_i , v)
            push!(buffer.GcInt2_i, ci_new)
        end            
    end
    temp1     = Gint[i]
    temp2     = GcInt2[i]
    Gint[i]   = buffer.Gint_i
    GcInt2[i] = buffer.GcInt2_i
    buffer.Gint_i   = temp1
    buffer.GcInt2_i = temp2
    return nothing
end


             
             
"""
    rk = upperTrapezoidal!(Gint, eInt, GcInt2, vActive, ieq)
    
Transform equations ieq..end to upper trapezoidal form and return the 
rank `rk` (`ieq <= rk <= length(Gint)`).
"""
function upperTrapezoidal!(Gint, eInt::Vector{Int}, GcInt, GcInt2, vActive::Vector{Bool}, ieq::Int)::Int
    j        = 0
    e        = 0
    vj       = 0
    pivot    = 0
    oldPivot = 1
    neq      = length(Gint)
    buffer   = Buffer()
    for k = ieq:neq
        # Search for a pivot in equations k:neq
        j = 0
        e = 0
        
        # Inspect only equations with one variable
        for k2 = k:neq
            if length(Gint[k2]) == 1
                j = getFirstActiveIndex(Gint, vActive, k2) 
                if j > 0
                    e = k2
                    break
                end
            end
        end
        
        if j == 0
            # Inspect only equations with two variables
            for k2 = k:neq
                if length(Gint[k2]) == 2
                    j = getFirstActiveIndex(Gint, vActive, k2) 
                    if j > 0
                        e = k2
                        break
                    end
                end
            end
        end
        
        if j == 0
            # Inspect all equations        
            for k2 = k:neq
                j = getFirstActiveIndex(Gint, vActive, k2) 
                if j > 0
                    e = k2
                    break
                end
            end
        end

        if j > 0
            # Extract infos
            pivot = GcInt2[e][j]
            vj    = Gint[e][j]
            deleteat!(Gint[e] , j)  # Remove variable vj in equation e (is added later to the front)
            deleteat!(GcInt2[e], j)
                    
            # Swap equations k and e
            if e != k
                swapEquations!(Gint, eInt, GcInt, GcInt2, k, e)
            end
        else
            # no pivot found: equations k:neq have no active variables (rank = k-1)
            return k-1
        end
        
        # Eliminate variable vj from equations k+1:neq
        for i = k+1:neq
            eliminateVariable!(Gint, GcInt2, k, i, pivot, vj, oldPivot, buffer)
        end
        oldPivot = pivot

        # Add variable vj in equation k at the front
        pushfirst!(Gint[k] , vj)
        pushfirst!(GcInt2[k], pivot)
        
        #println("... k = $k")
        #println("        Gint  = ", Gint)
        #println("        GcInt2 = ", GcInt2)
    end
    return neq
end



"""
    AvarRev = revertAssociation(Avar::Vector{Int})

Revert the association Vector `Avar[i] = j`, such that
`AvarRev[j] = i` (`Avar[i]=0` is allowed and is ignored).
"""
function revertAssociation(Avar::Vector{Int})::Vector{Int}
    AvarRev  = fill(0, length(Avar))
    for (i, value) in enumerate(Avar)
        if value != 0
            AvarRev[ value ] = i
        end
    end
    return AvarRev
end


needsVariableElimination(vProperty, eq_k::Vector{Int}) = findfirst(v->isEliminated(vProperty,v), eq_k) != nothing       


"""
    equationRemoved = simplifyOneEquation!(eq_k, eqc_k, AvarRev, vEliminated, vProperty)
    
Simplify one equation has much as possible. Return true, if equation is removed, otherwise return false.
"""
function simplifyOneEquation!(eq_k, eqc_k, AvarRev, vEliminated, vProperty)::Bool   
    while length(eq_k) > 1 && needsVariableElimination(vProperty, eq_k[2:end])            
        # Equation has more as one variable and at least one of the variables is zero, alias or negative alias variable    
        for j = 2:length(eq_k)
            v_j = eq_k[j]
            
            # Eliminate variable if possible
            if isNotEliminated(vProperty, v_j)
                continue
                
            elseif isZero(vProperty, v_j)
                deleteat!(eq_k , j)
                deleteat!(eqc_k, j)
                break
                
            else  # alias or negAlias
                if isAlias(vProperty, v_j)
                    v_add  = alias(vProperty, v_j)
                    vc_add = eqc_k[j]
                else
                    v_add  = negAlias(vProperty, v_j)
                    vc_add = -eqc_k[j]
                end
                
                # Check whether v_add appears in equation
                isPresent = false
                for i = 2:length(eq_k)
                    if i != j && eq_k[i] == v_add
                        isPresent = true
                        eqc_k[i] += vc_add
                        if eqc_k[i] == 0
                            if i < j 
                                deleteat!(eq_k , [i,j])
                                deleteat!(eqc_k, [i,j])
                            else
                                deleteat!(eq_k , [j,i])
                                deleteat!(eqc_k, [j,i])
                            end
                        else                                
                            deleteat!(eq_k , j)
                            deleteat!(eqc_k, j)
                        end
                        break
                    end
                end
                if isPresent
                    break
                else
                    eq_k[j]  = v_add
                    eqc_k[j] = vc_add
                end
            end                     
        end
    end 
    
         
    # Check if equation can be removed
    vk = eq_k[1] 
    if AvarRev[vk] == 0
        # vk is not a derivative of a variable -> it can be removed
        if length(eq_k) == 1  # equation k is a function of one variable
            # Variable is zero -> remove equation
            push!(vEliminated, vk)
            vProperty[vk] = 0
            empty!(eq_k)
            empty!(eqc_k)
            return true
            
        elseif length(eq_k) == 2 && abs(eqc_k[1]) == abs(eqc_k[2]) # equation k is a function of alias variables                               
            if eqc_k[1] > 0 && eqc_k[2] < 0  ||
               eqc_k[1] < 0 && eqc_k[2] > 0
                # Alias variable -> remove equation
                push!(vEliminated, vk)
                vProperty[vk] = eq_k[2]
                empty!(eq_k)
                empty!(eqc_k)                  
            else
                # Negative alias variable -> remove equation
                push!(vEliminated, vk)
                vProperty[vk] = -eq_k[2]
                empty!(eq_k)
                empty!(eqc_k)   
            end
            return true            
        end
    end   
    return false
end

        
function printLinearIntegerEquations(Gint, eInt, GcInt, var_name::Function; rk=(0,0,0))::Int
    ne = 0
    for i = 1:length(Gint)
        e = Gint[i]

        if length(e) > 0
            ne += 1
            # Construct equation
            vc     = GcInt[i][1]
            prefix = (vc == 1 ? "" : (vc == -1 ? "-" : string(vc) ))
            str    = "0 = " * prefix * var_name(e[1])
            for j = 2:length(e)
                v  = e[j]
                vc = GcInt[i][j]
                str = str *  (vc > 0 ? " + " : " - ")
                
                if abs(vc) != 1
                    str = str * string(abs(vc)) * "*"
                end
                str = str * var_name(v)
            end
            println(lpad(string(eInt[i]), 8), ": ", str)
        end
        if i == rk[1]
            println("    ---------- rk1")
        elseif i == rk[2]
            println("    ---------- rk2")
        elseif i == rk[3]
            println("    ---------- rk3")
        end
    end
    return ne
end
    
    
"""
    (vEliminated, vProperty, nvArbitrary, redundantEquations) = 
        simplifyLinearIntegerEquations!(G, eInt, GcInt, Avar)
        
Remove singularities of the **linear Integer equations** of a DAE system and simplify these equations as much as possible. 

The following **singularities** are fixed by this function:

- Equations that are redundant are removed.

- State constraint that are not structurally visible are transformed
  to a structurally visible form so that structural index reduction algorithms,
  such as the Pantelides algorithm, can handle these state constraints.
  
- Variables that can have an arbitrary value and do not appear in the remaining set of equations
  (so must be solved from the linear Integer equations) are set to zero.
  The calling function should print a warning message in such a case (if `nvArbitrary > 0`).


The following **simplifications** are performed recursively (`c` is an arbitrary Integer literal):

- An equation `c*v1 = 0` is removed and `v1` is replaced by zero in all occurrences
  of the linear Integer equations and all expressions with `v1` are simplified.
  
- An equation `c*v2 + c*v3 = 0` or `c*v2 - c*v3 = 0` is removed and
  `v2` is replaced by `v3` or by `-v3` in all occurrences of the linear Integer equations
  and all expressions with `v2` and `v3` are simplified.
  
  
Note:

- Input arguments `G, eInt, GcInt` are changed by a call of this function and
  represent the transformed linear Integer equations.
  The Abstract Syntax Tree (AST) of all these equations must be 
  replaced by new AST equations defined by the returned arguments.
  
- The number of operations in the transformed linear Integer equations is 
  guaranteed to be not larger as the original equations - with exception
  of the structurally visible state constraints, where the number of operations
  could be larger.
  
- Potential states (variables appearing differentiated) and derivatives
  of potential states are **not** eliminated by this function.


# Input arguments

- `G`: Bi-partite graph/incidence matrix of all equations. Typically: `G::Vector{Vector{Int}}`
       On entry, `G` is the original graph. On exit, the linear Integer equations of `G` are
       typically changed.
     
- `eInt::Vector{Int}`: `G[eInt]` are the linear Integer equations of `G`.
       On exit, `eInt` is reordered. If `length(G[eInt[i]]) = 0`, then the corresponding equation is eliminated.
     
- `GcInt`: `GcInt[i]` is the vector of Integer coefficients that are associated with the variables 
           of `G[eInt[i]]`. Typically: `GcInt::Vector{Vector{Int}}`.
           On exit, `GcInt` is reordered according to `eInt` and typically most of the coefficients have been changed.

- `Avar::Vector{Int}`: Defines the derivatives of the variables:
           `A[i] = if der(v_i) == v_k then k else 0`. This vector is not changed by `simplifyLinearIntegerEquations!`.


# Output arguments

- `vEliminated::Vector{Int}`: Variables that are eliminated.

- `vProperty::Vector{Int}`:  Defines the properties of the eliminated variables.
  These properties can be inquired with the following exported functions:
  
  o `isNotEliminated(vProperty, v)` - if variable v is not eliminated.
  o `isEliminated(vProperty, v)` - if variable v is eliminated.
  o `isZero(vProperty,v)` - if eliminated variable v is zero.
  o `isAlias(vProperty,v)` - if eliminated variable v is an alias variable v = v_alias
  o `isNegAlias(vProperty,v)` - if eliminated variable v is a negative alias variable v = -v_alias
  o `alias(vProperty,v)` - alias variable v_alias of eliminated variable v (v = v_alias).
  o `negAlias(vProperty,v)` - negated alias variable v_alias of eliminated variable v (v = -v_alias).
   
- `nvArbitrary::Int`: Variables `vEliminated[1:nvArbitrary]` are variables that can be arbitrarily set and that
              have been set to zero.
              
- `redundantEquations::Vector{Int}`: G[redundantEquations] are redundant equations that have been removed.


# Algorithm

The algorithm to remove the singularities is sketched in the paper:

- Otter, Elmqvist (2017):
  [Transformation of Differential Algebraic Array Equations to Index One Form](http://www.ep.liu.se/ecp/132/064/ecp17132565.pdf),
  section 5. Modelica'2017 Conference.
      
An error in this algorithm was fixed, the algorithm was improved to handle large equation systems
and to simplify equations as much as possible.      

# Main developer

[Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/), 
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)

"""
function simplifyLinearIntegerEquations!(G, eInt::Vector{Int}, GcInt, Avar::Vector{Int}; log::Bool=false, var_name::Function = v -> "???")
    nv = length(Avar)
    
    # Revert derivative association vector Avar
    AvarRev = revertAssociation(Avar)
    
    # Construct vActive1 and vActive2 for first and second transformation to upper trapezoidal form
    # (vActiveX[v] = false, if variable shall be ignored when transforming to upper trapezoidal form)
    #     vActive1[v] = true, if variable is not a potential state (does not appear differentiated)
    #                         and variable v must be solved from the linear Integer equations
    #     vActive2[v] = true, if variable is not a potential state (does not appear differentiated)
    vActive2 = [v==0 for v in Avar]
    vActive1 = copy(vActive2)
    for e in setdiff(collect(1:length(G)), eInt) # Equations that are not linear Integer equations
        for v in G[e]
            vActive1[v] = false
        end
    end
    
    # Save all variables with vActive1[v] = true in vShouldBeSolved
    vShouldBeSolved = findall(vActive1)
        

    # Construct the bi-partite graph of the linear Integer equations
    # (a deepcopy of the relevant part of G) 
    Gint = deepcopy(G[eInt])
    
    if log
        println("\n+++ Remove singularities")
        println("    Linear Integer equations:")
        printLinearIntegerEquations(Gint, eInt, GcInt, var_name)
        
        println("    Unknown variables:")
        unknowns = collect(DataStructures.OrderedSet([v for e in Gint for v in e]))
        for v in unknowns
            print(lpad(string(v),8), ": ", var_name(v))
            if v in vShouldBeSolved
                print("    (to be solved by equations)\n")
            elseif Avar[v] > 0
                print("    (potential state)\n")
            else 
                print("\n")
            end
        end  
    end
    
    
    # Construct a deepcopy of GcInt (this copy is modified by the transformation to upper trapezoidal form)
    GcInt2 = deepcopy(GcInt)
   
    # First transformation to upper trapezoidal form:
    #    Diagonal entries: Variables v that must be solved from the linear Integer equations 
    #println("\n... before upperTrapezoidal!:")
    #println("     Gint  = $Gint")
    #println("     GcInt2 = $GcInt2")     
    rk1 = upperTrapezoidal!(Gint, eInt, GcInt, GcInt2, vActive1, 1)
    
        # Eliminate variables that must be solved from the linear Integer equations
        # but are not diagonal entries (= variables can be arbitrarily set)
        vSolved     = [Gint[i][1] for i = 1:rk1]
        vEliminated = setdiff(vShouldBeSolved, vSolved)
        vProperty   = fill(IS_PRESENT, nv)    
        for v in vEliminated
            vProperty[v] = 0
        end
        nvArbitrary = length(vEliminated)
        
        if log
            println("\n    After first transformation to trapezoidal form (eliminate variables that must be solved):")
            printLinearIntegerEquations(Gint, eInt, GcInt2, var_name, rk=(rk1,0,0))
        end
       
    # Second transformation to upper trapezoidal form: Ignore potential states 
    rk2 = upperTrapezoidal!(Gint, eInt, GcInt, GcInt2, vActive2, rk1+1)
    
    if log
        println("\n    After second transformation to trapezoidal form (ignore potential states):")
        printLinearIntegerEquations(Gint, eInt, GcInt2, var_name, rk=(rk1,rk2,0))
    end    
    
    # Third transformation to upper trapezoidal form: All remaining variables are potential states
    fill!(vActive2, true)  
    rk3 = upperTrapezoidal!(Gint, eInt, GcInt, GcInt2, vActive2, rk2+1)
    #println("\n... after upperTrapezoidal!:")
    #println("     Gint  = $Gint")
    #println("     GcInt2 = $GcInt2")    
    #println("     rk1 = $rk1, rk2 = $rk2, rk3 = $rk3")
     
    if log
        println("\n    After third transformation to trapezoidal form (eliminate potential states):")
        printLinearIntegerEquations(Gint, eInt, GcInt2, var_name, rk=(rk1, rk2, rk3))
    end 
    
    # Simplify equations from equation rk2 upto equation 1
    for k = rk2:-1:1
        simplifyOneEquation!(Gint[k], GcInt2[k], AvarRev, vEliminated, vProperty)        
    end
 
    if log
        println("\n    After alias elimination:")
        printLinearIntegerEquations(Gint, eInt, GcInt2, var_name, rk=(rk1, rk2, rk3))
    end 
    
    #println("\n... after equation simplification:")
    #println("     Gint  = $Gint")
    #println("     GcInt2 = $GcInt2")    
    #println("     vEliminated = ", vEliminated)
    #println("     vProperty[vEliminated] = ",  vProperty[vEliminated])   

    # Update GcInt (use GcInt2[i] if it has not more unknowns as the corresponding GcInt[i] equation)
    equationsRemoved = false
    for i = 1:rk2
        if length(GcInt2[i]) <= length(GcInt[i])
            # Use transformed equation
            GcInt[i] = GcInt2[i]
        else
            # Use original equation
            Gint[i] = G[ eInt[i] ]
            equationRemoved = simplifyOneEquation!(Gint[i], GcInt[i], AvarRev, vEliminated, vProperty)
            equationsRemoved = equationsRemoved || equationRemoved
        end
    end
  
    if equationsRemoved
        # Simplify equations, until no equation is removed anymore
        equationRemoved = false
        while true
            for i = 1:rk2
                if length(Gint[i]) > 0
                    equationRemoved = simplifyOneEquation!(Gint[i], GcInt[i], AvarRev, vEliminated, vProperty)
                    if equationRemoved
                        break
                    end
                end
            end
            if !equationRemoved
                break
            end
        end
    end  
    
    # For constraint equations and for removed equations, use transformed equations
    for i = rk2+1:length(eInt)
        GcInt[i] = GcInt2[i]
    end
    
    if log
        println("\n    Final, simplified equations:")
        printLinearIntegerEquations(Gint, eInt, GcInt2, var_name, rk=(rk1, rk2, rk3))
    end 
    
    
    # Update G
    for (i,e) in enumerate(eInt)
        G[e] = Gint[i]
    end

    redundantEquations = rk3 < length(eInt) ? eInt[rk3+1:end] : Int[]
        
    return (vEliminated, vProperty, nvArbitrary, redundantEquations)
end



"""
    printSimplifiedLinearIntegerEquations(G, eInt, GcInt, vEliminated, vProperty, 
        nvArbitrary, redundantEquations, var_name::Function; printTest=false)
        
Print result of [`simplifyLinearIntegerEquations!`](@ref). 

Function `var_name(v)` returns the name of variable `v` as String.

If `printTest=true`, statements are printed that can be included in a Testset.
"""
function printSimplifiedLinearIntegerEquations(G, eInt, GcInt, vEliminated, vProperty, nvArbitrary, redundantEquations, var_name::Function; printTest=false)::Nothing    
    if nvArbitrary > 0
        println("\n    Variables that can be arbitrarily set and have been set to zero:")
        for v in vEliminated[1:nvArbitrary]
            println(lpad(string(v), 8), ": ", var_name(v), " = 0")
        end
    end
    
    
    if length(vEliminated) - nvArbitrary > 0
        println("\n    Variables that have been eliminated:")
        for v in vEliminated[nvArbitrary+1:end]
            print(lpad(string(v), 8), ": ", var_name(v), " = ")
            if isZero(vProperty, v)
                println("0")
            elseif isAlias(vProperty, v)
                println(var_name(alias(vProperty,v)))
            else
                println("-", var_name(negAlias(vProperty,v)))
            end
        end
    end
    
    
    if length(redundantEquations) > 0
        println("\n    Redundant equations that have been removed:")
        for e in redundantEquations
            println("      ", e)
        end
    end
    
    
    println("\n    Remaining transformed linear Integer equations:")
    ne = printLinearIntegerEquations(G[eInt], eInt, GcInt, var_name)
    if ne == 0
        println("      none (all linear Integer equations are removed)")
    end
    println()
    
    if printTest
        println("@test nvArbitrary == ", nvArbitrary)
        println("@test vEliminated == ", vEliminated)
        println("@test vProperty[vEliminated] == ", vProperty[vEliminated])
        println("@test redundantEquations == ", redundantEquations)           
        println("@test eInt  == ", eInt)
        println("@test G[eInt] == ", G[eInt])        
        println("@test GcInt == ", GcInt)
    end
    
    return nothing
end
