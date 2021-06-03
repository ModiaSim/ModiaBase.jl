# License for this file: MIT (expat)
# Copyright 2020-2021, DLR Institute of System Dynamics and Control
# Author: Martin Otter, DLR-SR
#
# Provide information about the structure of the ODE or DAE equations.

import DataFrames
import DataStructures
using  LinearAlgebra

export StateCategory, XD, XALG, XLAMBDA, XMUE
export ResidualCategory, FD, FC_ALG, FC_LOW_HIGH, FC_MUE

export LinearEquations, LinearEquationsIterator

export EquationInfo, StateElementInfo, update_equationInfo!, get_stateNames, get_x_table

export EquationInfoStatus, MANUAL, CODE_GENERATION, SOLVER_MODEL




"""
    @enum StateCategory

Type of a state.

| value     | description                                                                       |
|:----------|:----------------------------------------------------------------------------------|
| `XD`      | Differential variable (derivative of variable appears in equations)               |
| `XALG`    | Algebraic variable that is included in x                                          |
| `XLAMBDA` | Algebraic variable that is included in der_x                                      |
| `XMUE`    | Algebraic variable used for stabilization and included in der_x (exact value = 0) |
"""
@enum StateCategory XD=1 XALG XLAMBDA XMUE


"""
    @enum ResidualCategory

Type of a residual (see also [`StateCategory`](@ref) for the type of a state).

| value         | description                                                              |
|:--------------|:-------------------------------------------------------------------------|
| `FD`          | Differential equation (depends on XD, XLAMBDA, XMUE variables)           |
| `FC_ALG`      | Constraint equation that is purely algebraic (depends on XALG variables) |
| `FC_LOW_HIGH` | Constraint equation on lowest and on highest derivative level            |
| `FC_MUE`      | Constraint equation that is associated with a stabilizer XMUE            |
"""
@enum ResidualCategory FD=1 FC_ALG FC_LOW_HIGH FC_MUE



const niter_max = 20  # Maximum number of fixed-point iterations to solve A*x = b

"""
    leq = LinearEquations{FloatType}(vTear_names::Vector{String},  vTear_lengths::Vector{Int},
                                     nResiduals::Int, A_is_constant::Bool)

Define linear equation system "A*x=b" with `length(x) = sum(vTear_lengths)`.
If `A_is_constant = true` then `A` is a matrix that is constant after
initialization.

For details of its usage for code generation see [`LinearEquationsIterator`](@ref).
"""
mutable struct LinearEquations{FloatType <: Real}
    A_is_constant::Bool             # = true, if A-matrix is constant
    vTear_names::Vector{String}     # The names of the tearing variables (for error message)
    vTear_lengths::Vector{Int}      # Lengths of the vTear_names elements (length(vTear_value) = sum(vTear_lengths))    
    vTear_value::Vector{FloatType}  # Values of vTear variables (needed during execution of generated code)
    nResiduals::Int                 # Number of residuals
    
    A::Matrix{FloatType}
    b::Vector{FloatType}
    residuals::Vector{FloatType}    # Values of the residuals FloatType vector; length(residuals) = sum(residuals_length) = sum(vTear_lengths)   
    residual_value::AbstractVector  # Values of the residual variables ::Vector{Any}, length(residual_values) = nResiduals
    
    # Iteration status of for-loop
    mode::Int       # Operational mode (see function iterate below)
    niter::Int      # Current number of iterations in the fix point iteration scheme
    niter_max::Int  # Maximum number of iterations
    success::Bool   # = true, if solution of A*x = b is successfully computed
                    # = false, if solution is not computed; continue with fix point iteration
                    
    # For warning message if niter > niter_max
    inconsistentPositive::Vector{String}
    inconsistentNegative::Vector{String}   
    
    # Constructed during initialization
    residual_unitRanges::Vector{UnitRange{Int}} # residuals[residual_unitRanges[i]] = residual_value[i], if residual is a vector
    residual_indices::Vector{Int}               # residuals[residual_indices[i]] = residual_value[i], if residual is a scalar
    luA::LU{FloatType,Array{FloatType,2}}       # lu-Decomposition of A
    
    function LinearEquations{FloatType}(vTear_names::Vector{String}, vTear_lengths::Vector{Int},
                                        nResiduals::Int, A_is_constant::Bool) where {FloatType <: Real}
        @assert(length(vTear_names) > 0)
        @assert(length(vTear_names) == length(vTear_lengths))
        nx = sum(vTear_lengths)
        @assert(nx > 0)
        state = -1
        new(A_is_constant, vTear_names, vTear_lengths, zeros(FloatType,nx), nResiduals, 
            zeros(FloatType,nx,nx), zeros(FloatType,nx), zeros(FloatType,nx),
            Vector{Any}(undef,nResiduals), -4, -1, niter_max, false, String[], String[],
            fill(0:0,nResiduals), fill(0,nResiduals))
    end
end
LinearEquations(args...) = LinearEquations{Float64}(args...)
    


"""
    leqIterator = LinearEquationsIterator{FloatType}(leq::LinearEquations)

Return instance of a struct to iterate over LinearEquations with Base.iterate(leq::LinearEquationsIterator)
with the following for-loop:

```julia
function getDerivatives!(_der_x, _x, _m, _time)::Nothing
    _leq::Union{Nothing,LinearEquations{FloatType}} = nothing
    ...
    _leq = _m.linearEquations[<nr>]   # leq::LinearEquations{FloatType}
    for _leq_evaluate = LinearEquationsIterator(_leq, _m.isInitial, _m.time) 
        if !_leq_evaluate
            break
        end
        v_tear1 = _leq.vTear_value[1:3]
        v_tear2 = _leq.vTear_value[4]
        ...
        v_solved1 = f(v_tear1, v_tear2, ..., positive(v_tear1,.., _leq))
        v_solved2 = f(v_tear1, v_tear2, ..., v_solved1)
        ...
        _leq.residuals[..] = < residual equation(v_tear1,...) > 
    end
    ...
end
```
"""
struct LinearEquationsIterator{FloatType <: Real}
    leq::LinearEquations{FloatType}
    isInitial::Bool
    time
end

const ncalls = [0]

function Base.iterate(iterator::LinearEquationsIterator{FloatType}, mode::Int = -3) where {FloatType <: Real}
    #=
        for item in iterator
           # body
        end
        
    is translated to
    
        next = iterate(iterator)
        while next != nothing
            (item, mode) = next
            # body
            next = iterate(iterator, mode)
        end
    
    
    iterate(iterator,mode=0) to solve "residuals = A*x - b" where
    A and b can be functions of positive(c_i'*x - d_i). In this case this system
    is solved by a fixed point iteration scheme, that is "A*x = b" is solved
    until potentially present positive(c_i'*x - d_i) calls are consistent to x.
    Iteration is terminated after 20 iterations if no convergence is found.
    In this case a warning message is triggered and simulation is continued
    (it might be that simulation is still successful, even if "x" and
    positive(..) are temporarily not consistent to each other).

    Note, the current values of A,x,b,residuals are stored in iterator.leq.
    If A is fixed after initialization (leq.A_is_constant = true), then A is computed
    only once at initialization, the LU decomposition of A is stored in leq
    and used in subsequent calls to solve the equation system.
    
    
        On input (leq = iterator.leq; x = leq.vTear_value):
        
        mode = -3: # Initialize fixed point iteration
                   leq.niter   = 0
                   leq.success = false 
                   goto mode = 0
                   
        mode =  0: # Start next (fixed point) iteration 
                   x .= 0
                   leq.mode = 0
                   return (true, 1)  # next mode = 1 (compute "b" after the next iteration)
                
        mode = 1: b    = -residuals  # Compute "b". 
                  x[1] = 1.0
                  leq.mode = 1
                  return (true, 2)   # next mode = 2 (compute A[:,1] after the next iteration)
                
        mode > 1: j      = mode - 1
                  A[:,j] = residuals + b  # Compute A[:,j]. 
                  x[j]   = 0.0
                  if j < length(x)
                      x[j+1]   = 1.0
                      leq.mode = mode
                      return (true, mode+1)  # next mode = mode + 1 (compute A[:,j+1] after the next iteration)
                  elseif j == length(x)
                      x = A\b            # Solve linear equation system A*x = b for x.
                      leq.mode    = -1   # Compute all variables IN the next iteration as function of x
                      leq.niter  += 1    # Increment number of iterations to solve A*x = b
                      leq.success = true # Terminate for-loop at the beginning of the next iteration,
                                         # provided no positive(..) call changes its value.
                                         # (leq.success is set to false in positive(..), if the return value changes).
                      return (true, -2)  # next mode = -2 (either terminate for-loop or re-initialize the next solver-iteration)
                  end
                  
       mode = -2: if leq.success
                     return (false, 0)   # Terminate for-loop
                  elseif leq.niter > 20
                     <warning>
                     return (false, 0)   # Terminate for-loop
                  end
                  goto mode = 0          # Re-initialize next solver-iteration
    =# 
    
    #=
    println("... 1, mode = $mode, time = $(iterator.time)")
    if iterator.isInitial
        ncalls[1] = 0
    else
        ncalls[1] += 1
        if ncalls[1] > 20000
        error("... too many calls")
        end
    end
    =#
    
    
    leq = iterator.leq
    if mode == -3 
        # Initialize fixed point iteration
        leq.niter   = 0
        leq.success = false 
        empty!(leq.inconsistentPositive)
        empty!(leq.inconsistentNegative)        
        mode = 0
        
    elseif mode == -2
        # Either terminate fixed point iteration, or start a new iteration
        if leq.success
            return (false,0)
        elseif leq.niter > niter_max
            str = ""
            if length(leq.inconsistentPositive) > 0
                str = str * "positive(expr) is inconsistent for expr = $(leq.inconsistentPositive)."
            end
            if length(leq.inconsistentNegative) > 0
                if length(str) > 0
                    str = str * "\n"
                end
                str = str * "negative(expr) is inconsistent for expr = $(leq.inconsistentNegative)."
            end
            @warn "At time = $(iterator.time), no consistent solution found for mixed linear equation system.\n" * 
                  "Simulation is continued although some variables might not be correct at this time instant.\n$str"
            return (false,0)
        end
        mode = 0
    end
    
    x = leq.vTear_value
    if mode == 0
        # Re-initialize iteration variables and compute b-vector after next iteration
        x .= 0
        leq.mode = 0
        return (true,1)
    end

    nx  = length(x)
    A   = leq.A
    b   = leq.b
    residuals           = leq.residuals
    residual_value      = leq.residual_value
    residual_unitRanges = leq.residual_unitRanges
    residual_indices    = leq.residual_indices

    if mode < 1 || mode > nx+1
        @goto ERROR
    end
    
    if iterator.isInitial && mode == 1
        # Construct unit ranges for the residual variables vector to copy values into the residuals vector
        j = 1
        for (i, res_value) in enumerate(residual_value)
            if typeof(res_value) <: Number
                residual_indices[i] = j
                j = j+1
            else
                len = length(res_value)
                if length(res_value) == 0
                    error("Residual variable value $i has length zero")
                end
                k = j+len-1
                residual_unitRanges[i] = j:k
                j = k+1
            end
        end
        if j-1 != nx
            k = j-1
            #@error "The length of the residuals vector (= $k) is not equal to the length of the vTear_value vector (= $nx)"
        end
    end

    # Copy residual variable values to residuals vector
    index = 0
    for (i, res_value) in enumerate(residual_value)
        index = residual_indices[i]
        if index > 0
            residuals[index] = res_value
        else
            residuals[residual_unitRanges[i]] = res_value
        end
    end

    if !leq.A_is_constant || iterator.isInitial  # A is not constant or A is constant and isInitial = true
        if mode == 1
            # Terminating code for mode = 0 (residuals = A*x - b -> b = -residuals)
            for i = 1:nx
                b[i] = -residuals[i]
            end
        else
            # Terminating code for mode = 1..nx (residuals = A*x - b -> A[:,j] = residuals + b)
            j = mode-1
            for i = 1:nx
                A[i,j] = residuals[i] + b[i]
            end
            x[j] = 0
        end
    
        if mode <= nx
            # Start code for mode = 1..nx
            x[mode] = 1
            leq.mode = mode
            return (true, mode+1)
        end
        
        # mode == nx+1; Solve linear equation system
        if length(x) == 1       
            x[1] = b[1]/A[1,1]
            if !isfinite(x[1])
                error("Linear scalar equation system is singular resulting in: ", leq.vTear_names[1], " = ", x[1])
            end
        else
            x .= b
            leq.luA = lu!(A)
            ldiv!(leq.luA, x)
        end
        
    elseif leq.A_is_constant && !iterator.isInitial && mode == 1  # isInitial=false, LU decomposition of A is available in leq.luA    
        for i = 1:nx
            x[i] = -residuals[i]
        end 
            
        # Solve linear equation system
        if length(x) == 1
            x[1] = x[1]/A[1,1]
            if !isfinite(x[1])
                error("Linear scalar equation system is singular resulting in: ", leq.vTear_names[1], " = ", x[1])
            end
        else        
            ldiv!(leq.luA, x)
        end
        
    else
        @goto ERROR
    end   
   
    # Linear equation system solved
    leq.mode    = -1   # Compute all variables IN the next iteration as function of x
    leq.niter  += 1    # Increment number of iterations to solve A*x = b
    leq.success = true # Terminate for-loop at the beginning of the next iteration,
                       # provided no positive(..) call changes its value.
                       # (leq.success is set to false in positive(..), if the return value changes).
    return (true, -2)  # next mode = -2 (either terminate for-loop or re-initialize the next solver-iteration)            

    @label ERROR
    @error "Should not occur (Bug in file ModiaBase/src/EquationAndStateInfo.jl,\nmode = $mode; vTear = $(leq.vTear_names))."
end




"""
    @enum EquationInfoStatus

Status of an EquationInfo instance:

- `MANUAL`: Is defined manually. The following variables in x_info::Vector{StateElementInfo}
   are **not** defined:
   x_names_julia, der_x_names_julia, length, unit, startIndex.
   Also variables nx and x_infoByIndex are not defined.
   With function [`update_equationInfo!`](@ref), the variables length, unit, startIndex
   in x_info::Vector{StateElementInfo} are computed, as well as nx and x_infoByIndex.
   
- `CODE_GENERATION`: Is constructed during code generation with getSortedAndSolvedAST(..).
  The following variables in x_info::Vector{StateElementInfo}
  are **not** defined: startIndex.
  Also variables nx and x_infoByIndex are not defined.  
  With function [`update_equationInfo!`](@ref), the variables startIndex
  in x_info::Vector{StateElementInfo} are computed, as well as nx and x_infoByIndex.
  
- `SOLVER_MODEL`: Is used during simulation in a SolverModel. With function
  [`update_equationInfo!`](@ref), missing variables are constructed depending
  on the information given with `MANUAL` or `CODE_GENERATION` and the actual
  modelValues instance (here unit and length information is available). 
"""
@enum EquationInfoStatus MANUAL=1 CODE_GENERATION SOLVER_MODEL



"""
    xe_info = StateElementInfo(...)

Return an instance of the mutable struct `StateElementInfo` that defines the information
for one element of the state vector. There are three constructors:

- Default constructor (all variables under section Arguments are given;
  used to write/read the complete information).
  
- All variables under section Arguments are given with exception of startIndex
  (used for `EquationInfoStatus = CODE_GENERATION`).

- StateElementInfo(x_name, der_x_name, stateCategory; fixed=nothing, nominal=NaN, unbounded=false)
  (used for `EquationInfoStatus = MANUAL`):

# Arguments
- x_name: Name of x-element or "" if no name (if stateCatebory = XLAMBDA or XMUE)
- x_name_julia: Julia name of x-element in getDerivatives!/getResiduals! function
  or "" if not needed (since no code generation).
- der_x_name: Name of der_x-element or "" if either `der(x_name)` or if no name
  (if stateCategory = XALG).
- der_x_name_julia: Julia name of der_x-element in getDerivatives!/getResiduals! function
  or "" if not needed (since no code generation).
- stateCategory::StateCategory: Category of the state
- length: length of x-element (or -1 if not yet known) 
- unit: unit of XD, XALG (x_name) or XLAMBDA, XMUE (der_x_name) or "" if not yet known)
- fixed: false (= guess value) or true (= not changed by initialization).
         Only relevant for ode=false, otherwise ignored.
- nominal: Nominal value (NaN if determined via start value)
- unbounded: false or true
- startIndex: start index of x-element with respect to x-vector 
              or -1 if not yet known.
"""
mutable struct StateElementInfo
    x_name::String                # Modia name of x-element or "" if no name (for )
    x_name_julia                  # Julia name of x-element in getDerivatives! function
                                  # or "", if not needed (since no code generation).
    der_x_name::String            # Modia name of der_x-element or "" if either "der(x_name)" or if no name,
    der_x_name_julia              # Julia name of der_x-element in getDerivatives! function
                                  # or not needed (since no code generation)
    stateCategory::StateCategory  # category of the state
    length::Int                   # length of x-element (or -1 if not yet known)       
    unit::String                  # unit of x-element as string (or "" if not yet known)
    fixed::Bool                   # false or true
    nominal::Float64              # nominal value (NaN if determined via start value)
    unbounded::Bool               # false or true
    startIndex::Int               # start index of x-element with respect to x-vector 
                                  # or -1 if not yet known.
end

# Constructor for code-generation 
StateElementInfo(x_name, x_name_julia, der_x_name, der_x_name_julia,
                 stateCategory, length, unit, fixed, nominal, unbounded) = StateElementInfo(
                 x_name, x_name_julia, der_x_name, der_x_name_julia,
                 stateCategory, length, unit, fixed, nominal, unbounded, -1)

# Constructor for reading StateElementInfo after code-generation (x_name_julia and der_x_name_julia are not included
StateElementInfo(x_name, der_x_name, 
                 stateCategory, length, unit, fixed, nominal, unbounded) = StateElementInfo(
                 x_name, :(), der_x_name, :(),
                 stateCategory, length, unit, fixed, nominal, unbounded, -1)

# Constructor for manually generated equation info.					 
StateElementInfo(x_name, der_x_name, stateCategory=XD; fixed=false, nominal=NaN, unbounded=false) =
    StateElementInfo(x_name, :(), der_x_name, :(), stateCategory, -1, "", fixed, nominal, unbounded, -1)


function Base.show(io::IO, xe_info::StateElementInfo)
    print(io, "ModiaBase.StateElementInfo(")
    show( io, xe_info.x_name)
	#print(io, ",")		
    #show( io, xe_info.x_name_julia)
	print(io, ",")			
    show( io, xe_info.der_x_name)
	#print(io, ",")			
    #show( io, xe_info.der_x_name_julia)
    print(io, ",ModiaBase.", xe_info.stateCategory)
    print(io, ",", xe_info.length)
	print(io, ",")	
    show( io, xe_info.unit)
    print(io, ",", xe_info.fixed)
    print(io, ",", xe_info.nominal)
    print(io, ",", xe_info.unbounded)
    #if xe_info.startIndex > 0
    #    print(io, ",", xe_info.startIndex)
    #end
    print(io, ")")
    return nothing
end


"""
    x_table = get_x_table(x_info::Vector{StateElementInfo}
    
Return the state element info as DataFrames.DataFrame table, for example
to print it as:

```
show(x_table, allrows=true, allcols=true, summary=false, eltypes=false)
```
"""
function get_x_table(x_info::Vector{StateElementInfo})
    x_table = DataFrames.DataFrame(name=String[], length=Int[], unit=String[], fixed=Bool[], nominal=Float64[], unbounded=Bool[])
    
    for xe_info in x_info
        push!(x_table, (xe_info.x_name, xe_info.length, xe_info.unit, xe_info.fixed, xe_info.nominal, xe_info.unbounded))
    end
    
    return x_table
end


"""
    eqInfo = EquationInfo(;
                status               = MANUAL,
                ode                  = true, 
                nz                   = 0,
                x_info               = StateElementInfo[],
                residualCategories   = ResidualCategory[],
                linearEquations      = Tuple{Vector{String},Bool}[],             
                vSolvedWithFixedTrue = String[],
                defaultParameterAndStartValues = nothing,
                ResultType = nothing,
                ResultTypeHasFloatType = false)

Return instance `eqInfo` that defines the information for the equation system.

# Arguments
- status::EquationInfoStatus: Defines the variables that have a value.
- ode: = true if ODE-interface (`getDerivatives!`), 
       = false if DAE-interface (`getResiduals!`).
- nz: Number of zero crossing functions.
- x_info: Vector of StateElementInfo elements provding info for every x-element     
- residualCategories: If ode=true, length(residualCategories) = 0.
             If ode=false: residualCategories[i] is the `ResidualCategory`](@ref) of x-element "i".
- `linearEquations::Vector{Tuple{Vector{String},Vector{Int},Int,Bool}}`: 
               linearEquations[i] defines a 
               ModiaBase.LinearEquations system, where the first tuple value 
               is a vector of the names of the unknowns, the second tuple value
               is a vector with the lengths of the unknowns, the third tuple value is the number
               of residuals and the fourth tuple value 
               defines whether the coefficient matrix A
               has only constant entries (=true) or not (=false).             
- `vSolvedWithFixedTrue::Vector{String}`: Vector of variables that are computed
               from other variables and have `fixed=true`. During initialization
               it is checked whether the calculated values and the start values of
               these variables agree. If this is not the case, an error is triggered.   
- `defaultParameterAndStartValues`: Dictionary of default paramter and default start values.
- `ResultType::AbstractModelValues`: ModelValues type that shall be used for the result generation
               (= ModelValues struct without parameters).  
- `ResultTypeHasFloatType::Bool`: = false, if `ResultType` is not parameterized.
                                  = true, if `ResultType` has parameter `FloatType`.               
"""
mutable struct EquationInfo
    status::EquationInfoStatus
    ode::Bool                                      # = true if ODE model interface, otherwise DAE model interface
    nz::Int                                        # Number of crossing functions
    x_info::Vector{StateElementInfo}
    residualCategories::Vector{ResidualCategory}   # If ode=true, length(residualCategories) = 0
                                                   # If ode=false, residualCategories[j] is the ResidualCategory of residual[j] 
    linearEquations::Vector{Tuple{Vector{String},Vector{Int},Int,Bool}}
    vSolvedWithFixedTrue::Vector{String}
    nx::Int                                        # = length(x) or -1 if not yet known
    x_infoByIndex::Vector{Int}                     # i = x_infoByIndex[j] -> x_info[i] 
                                                   # or empty vector, if not yet known.  
    x_dict::DataStructures.OrderedDict{String,Int} # x_dict[x_name] returns the index of x_name with respect to x_info                                                  
    defaultParameterAndStartValues::Union{AbstractDict,Nothing}
    ResultType    
    ResultTypeHasFloatType::Bool
end
 
 
EquationInfo(; status                = MANUAL,
               ode                   = true, 
               nz                    = 0,
               x_info                = StateElementInfo[],
               residualCategories    = ResidualCategory[],
               linearEquations       = Tuple{Vector{String},Vector{Int},Int,Bool}[],
               vSolvedWithFixedTrue  = String[],
               nx                    = -1,
               x_infoByIndex         = Int[],
               defaultParameterAndStartValues::Union{AbstractDict,Nothing} = nothing,
               ResultType = nothing,
               ResultTypeHasFloatType = false) = 
               EquationInfo(status, ode, nz, x_info, 
                            residualCategories, linearEquations, 
                            vSolvedWithFixedTrue, nx, x_infoByIndex,
                            DataStructures.OrderedDict{String,Int}(),
                            defaultParameterAndStartValues,
                            ResultType, ResultTypeHasFloatType)


"""
    updateEquationInfo!(eqInfo::EquationInfo)
    
Set eqInfo.x_dict, eqInfo.nx and eqInfo.x_info[:].startIndex
"""
function updateEquationInfo!(eqInfo::EquationInfo)::Nothing
    x_dict = eqInfo.x_dict
    startIndex = 1
    for (i, xi_info) in enumerate(eqInfo.x_info)
        x_dict[xi_info.x_name] = i
        xi_info.startIndex = startIndex
        startIndex += xi_info.length
    end
    eqInfo.nx = startIndex - 1
    return nothing
end

get_x_startIndexAndLength(eqInfo::EquationInfo, name::String) = begin
    xe_info = eqInfo.x_info[ eqInfo.x_dict[name] ]
    (xe_info.startIndex, xe_info.length)
end
            
            
function Base.show(io::IO, eqInfo::EquationInfo; indent=4)
    indentation  = repeat(" ", indent)
    indentation2 = repeat(" ", 2*indent)
    indentation3 = repeat(" ", 3*indent)	
    println(io, "ModiaBase.EquationInfo(")
    println(io, indentation2, "status = ModiaBase.", eqInfo.status, ",")
    println(io, indentation2, "ode = ", eqInfo.ode, ",")
    if eqInfo.nz > 0
        println(io, indentation2, "nz = ", eqInfo.nz, ",")
    end    
    print(io, indentation2, "x_info = ModiaBase.StateElementInfo[")
    for (i, xe_info) in enumerate(eqInfo.x_info)
        if i == 1
            print(io, "\n", indentation3)
        else
            print(io, ",\n", indentation3)
        end
        show(io, xe_info)
    end
    print(io,"]")
    
    if length(eqInfo.residualCategories) > 0
        print(io, ",\n", indentation2, "residualCategories = [")
        for (i, rcat) in enumerate(eqInfo.residualCategories)
            if i > 1
                print(io, ",")
            end
            show(io, rcat)
        end
        print(io,"]")
    end

    leqs = eqInfo.linearEquations 
    if length(leqs) > 0
        println(io, ",\n", indentation2, "linearEquations = [")
        for (i, leq) in enumerate(leqs)
            print(io, indentation3, "(", leq[1], ",\n", 
                      indentation3, " ", leq[2], ",\n", 
                      indentation3, " ", leq[3], ", ", leq[4], ")")
            if i < length(leqs)
                print(io, ",\n")
            end
        end
        print(io, "]")
    end


    if length(eqInfo.vSolvedWithFixedTrue) > 0
        print(io, ",\n", indentation2, "vSolvedWithFixedTrue = ")
        show(io, eqInfo.vSolvedWithFixedTrue)
    end
    
    if eqInfo.nx > 0
        print(io, ",\n", indentation2, "nx = ", eqInfo.nx)
    end
    
    if length(eqInfo.x_infoByIndex) > 0
        print(io, ",\n", indentation2, "x_infoByIndex = ")
        show(io, eqInfo.x_infoByIndex) 
    end
    
    if !isnothing(eqInfo.defaultParameterAndStartValues)
        print(io, ",\n", indentation2, "defaultParameterAndStartValues = ")
        show(io, eqInfo.defaultParameterAndStartValues, indent=12, finalLineBreak=false) 
    end
    
    if !isnothing(eqInfo.ResultType)
        print(io, ",\n", indentation2, "ResultType = ")
        show(io, eqInfo.ResultType) 
        
        if eqInfo.ResultTypeHasFloatType
            print(io, ",\n", indentation2, "ResultTypeHasFloatType = ")
            show(io, eqInfo.ResultTypeHasFloatType) 
        end
    end    
    
    
    println(io, "\n", indentation, ")")
    return nothing
end



"""
    names = get_stateNames(equationInfo::EquationInfo)
    
Return the names of the states defined in `equationInfo` as a Vector of strings.
"""
get_stateNames(eqInfo::EquationInfo) = String[xi_info.x_name for xi_info in eqInfo.x_info]




#=
"""
    update_equationInfo!(eqInfo::EquationInfo)
    
Update equationInfo as needed for [`EquationInfoStatus`](@ref)` = SOLVER_MODEL`.
"""
function update_equationInfo!(eqInfo::EquationInfo, modelValues::ModiaBase.AbstractModelValues)::Nothing
    x_info = eqInfo.x_info
    
    if eqInfo.status == MANUAL
        # Determine length and unit
		count = 0
        for xi_info in x_info
            if xi_info.x_name != "" # XD or XALG
                (component,name) = ModiaBase.get_modelValuesAndName(modelValues, xi_info.x_name)
            elseif xi_info.der_x_name != ""  # XLAMBDA or XMUE
                (component,name) = ModiaBase.get_modelValuesAndName(modelValues, xi_info.der_x_name)
            else
				# Pure algebraic with dummy equation der(x) = -x, x(0) = 0
				xi_info.length = 1
				count = count+1
				continue
            end
            
	        xi_type        = fieldtype(typeof(component), name)
			xi_numberType  = ModiaBase.numberType(xi_type)			
            xi_info.unit   = replace(string(unit(xi_numberType)), " " => "*")
            xi_info.length = Base.length(xi_type)  # Base.length(::Type{<:Number})=1 is defined in ModiaBase, since not defined in Base
        end
		if count == 1 && length(x_info) > 1 || count >  1
			error("Error in update_equationInfo!: x_info is wrong.\n",
				  "x_info = $x_info")
		end
    end          
   
    # Set startIndex in x_info and compute nx
    startIndex = 1
    for xi_info in x_info
        xi_info.startIndex = startIndex
        startIndex += xi_info.length
    end
    nx = startIndex - 1
    eqInfo.nx = nx
            
    # Set x_infoByIndex::Vector{Int}
    eqInfo.x_infoByIndex = zeros(Int,nx)
    i1 = 0
    i2 = 0
    for (i,xi_info) in enumerate(x_info)
        i1 = xi_info.startIndex
        i2 = i1 + xi_info.length - 1
        for j = i1:i2
            eqInfo.x_infoByIndex[j] = i
        end
    end
    
    eqInfo.status = SOLVER_MODEL
		
    return nothing
end
=#

