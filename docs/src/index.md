# ModiaBase.jl Documentation

[ModiaBase](https://github.com/ModiaSim/ModiaBase.jl) provides functions to support the transformation of a
Differential Algebraic Equation system (DAE)

```math
0 = f_{dae}(\frac{dx_{dae}}{dt},x_{dae},w,t) \tag{1}
```

to an explicit Ordinary Differential Equation system (ODE)

```math
\begin{aligned}
  \frac{dx_{ode}}{dt} &= f_{x}(x_{ode},t) \\
         (w, x_{dae}) &= f_{w}(x_{ode},t)
\end{aligned} \tag{2}
```

where ``x=x(t), w=w(t)`` are vector valued functions of the scalar
variable ``t`` (= usually time). ``x_{dae}`` is the DAE state vector,
``x_{ode}`` is the ODE state vector - a subset of ``x_{dae}``, and
``w`` are purely algebraic variables that do not appear differentiated in the DAE.
The equations are hereby represented as a vector of Julia expressions,
that is of an Abstract Syntax Tree (AST).

These functions are used by package [Modia](https://github.com/ModiaSim/Modia.jl),
but can also be utilized in another context. Especially the following functionality is provided:

- Simplify linear Integer equations (many equations of object-oriented models are linear Integer equations and can be pre-processed exactly)
  - to remove alias variables and equations,
  - to remove redundant equations,
  - to provide definite values for variables that can have arbitrary values if this makes sense,
  - to make state constraints structurally visible.
  
- Find a variable assignment of an equation system, in order
  to transform the equation system in a directed graph that can be further
  processed.
  
- Find the strong components in a directed graph (with the algorithm of Tarjan)
  to determine algebraic equation systems that must be solved together.

- Sort an equation system (= transform to Block Lower Triangular form), 
  to determine the order in which the equations have to be evaluated.
  
- Reduce the dimension of algebraic equation systems by tearing.
 
- Find equations that need to be differentiated one or more times (with the algorithm of Pantelides)
  in order that the DAE can be transformed to an ODE.

- Analytically differentiate the found equations.

- Statically select ODE states and transform to ODE form
  (hereby identifying linear equation systems that must be solved during simulation).
  
Transformation from a DAE to an ODE form is (currently) performed if no nonlinear-algebraic equations
appear and the ODE-states can be statically selected.

Array variables and array equations are kept (they are not "flattened" in single elements).
However, DAE forms that require to differentiate array equations, are not yet supported.

The following extensions are planned (internal prototypes are available):

- Full support of array equations.
- If transformation to an ODE is not possible with the algorithms above,
  transformation to a special index 1 DAE, that
  can be simulated with standard DAE solvers (such as Sundials IDA).


## Installation

The package is registered and is installed with (Julia 1.7 is required):

```julia
julia> ]add ModiaBase
```


## Release Notes

### Version 0.9.0-dev

Non-backwards compatible improvements

- Parameter values in the code are now type cast to the type of the parameter value from the 
  `@instantiatedModel(..)` call. The benefit is that access of parameter values in the code is type stable
  and operations with the parameter value are more efficient and at run-time no memory is allocated.
  Existing models can no longer be simulated, if parameter values provided via `simulate!(.., merge=xx)` are not
  type compatible to their definition. For example, an error is thrown if the @instantedModel(..) uses a Float64 value and the
  `simulate!(.., merge=xx)` uses a `Measurement{Float64}` value for the same parameter.
  
Other improvements

- Hierarchical names in function calls supported (e.g. `a.b.c.fc(..)`). 

- Functions can return multiple values, e.g. `(tau1,tau2) = generalizedForces(derw1, derw2)`.

- Large speedup of symbolic transformation, if function depends on many input (and output) arguments 
  (includes new operator `implicitDependency(..)`).
  
- Support for StaticArrays variables (the StaticArrays feature is kept in the generated AST).
  
- Support for Array variables (especially of state and tearing variables)
  where the dimension can change after @instantiateModel(..)

- Included DAE-Mode in solution of linear equation system (if DAE integrator is used and all unknowns of a linear
  equation system are part of the DAE states, solve the linear equation system during continuous integration
  via DAE solver (= usually large simulation speed-up, for larger linear equation systems)  
 
Bug fixes

- Do no longer expand the unit macro in the AST, such as `u"N"`, because otherwise `logCode=true` results in wrong code
  (previously, a `u"N"`  definition in the model was displayed in the code as `N` which is not correct Julia code).


### Version 0.8.1

- Update Project.toml, Manifest.toml, README.md


### Version 0.8.0

- Require Julia 1.7
- Upgrade Manifest.toml to version 2.0
- Update Project.toml/Manifest.toml

### Version 0.7.8

- Tests of TestDifferentiate.jl corrected to comply with DiffRules > 1.0
- Scaling introduced to improve numerics when constructing A-matrix of linear equation system.


### Version 0.7.7

- Bug fixed when selecting RecursiveFactorization.jl


### Version 0.7.6

- Fixed bug in StateSelection.jl: If unitless=true, no unit is associated with the tearing variable.

- Solve linear equation systems optionally with [RecursiveFactorization.jl](https://github.com/YingboMa/RecursiveFactorization.jl) 
  instead of the default `lu!(..)` and `ldiv!(..)`.

- Project.toml: Changed DiffRules from "~1.0" to "1", since issue with "1.2.1" 
  (leading to an error in runtests) seems to be fixed.
  
- Project.toml: Added version 1 of MonteCarloMeasurements.
  
- Updated used packages.

- Tutorial slightly improved.


### Version 0.7.5

- Added a restriction, so that DiffRules 1.0.2 is used, instead of 1.2.1 (which leads to an error in the test suite).


### Version 0.7.4

- showCodeWithoutComments(code): Bug corrected to only remove comments and not other code
  (ModiaLang.@instantiateModel(..., logCode=true, ...) gave wrong output).
  
- Used packages updated


### Version 0.7.3

- Speed improvements for structural and symbolic algorithms.

- Added support for state events, time events and synchronous operators
  (positive(), Clock(), after(), pre(), previous(), hold(), initial(), terminal()) 

- Added support for mixed linear equation systems having Real and Boolean unknowns.

- Simplified code for linear equation systems (while-loop instead of for-loop).

- Added TimerOutputs @timeit instrumentation to the solution of linear equation systems.


### Version 0.7.2

- Support of parameters as hierarchical named tuples.

- Support of array comprehensions.

- Support of array end (e.g. A[3:end])

- If one equation cannot be solved for one unknown (e.g. since function call),
  try to solve it as linear equation system.
  
- If variables with init values are explicitly solved for, print warning message
  only if log = true (in TinyModia.simulate! an error occurs, if the init value
  cannot be respected).


### Version 0.7.1

- Due to version conflicts, added version 0.17 of DataStructures in compat.


### Version 0.7.0

- Initial version, based on code developed for Modia 0.6 and ModiaMath 0.6.


## Main developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)