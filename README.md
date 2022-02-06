# ModiaBase
 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://modiasim.github.io/ModiaBase.jl/stable/)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/ModiaSim/ModiaBase.jl/blob/master/LICENSE.md)

ModiaBase is part of [ModiaSim](https://modiasim.github.io/docs/). It is usually used via [Modia](https://github.com/ModiaSim/Modia.jl).
The [ModiaBase documentation](https://modiasim.github.io/ModiaBase.jl/stable/) provides details of the algorithms and how to use them.

ModiaBase provides basic algorithms and functionality that is needed for
equation-based modeling to transform a (potentially high-index) Differential-Algebraic Equation system (DAE),
to an Ordinary Differential Equation system in state space form (ODE).
It is used by [Modia](https://github.com/ModiaSim/Modia.jl),
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
  

## Installation
 
Typically, a user installs [Modia](https://github.com/ModiaSim/Modia.jl) and does not need
to install ModiaBase separately. If needed, ModiaBase is installed with (Julia 1.7 is required):

```julia
julia> ]add ModiaBase
```

## Main Developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).
  

License: MIT (expat)
