# ModiaBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://modiasim.github.io/ModiaBase.jl/stable/)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/ModiaSim/ModiaBase.jl/blob/master/LICENSE.md)

ModiaBase provides basic algorithms and functionality that is needed for
equation-based modeling to transform a (potentially high-index) Differential-Algebraic Equation system (DAE),
to an Ordinary Differential Equation system in state space form (ODE).
It is used by [TinyModia](https://github.com/ModiaSim/TinyModia.jl),
but can also be utilized in another context. Especially the following functionality is provided:

- Find a variable assignment of an equation system, in order
  to transform the equation system in a directed graph that can be further
  processed.

- Find the strong components in a directed graph (with the algorithm of Tarjan)
  in order to find algebraic equation systems that must be solved together.

- Sort an equation system (= transform to Block Lower Triangular form), in order
  determine the order in which the equations have to be evaluated.

- Simplify linear Integer equations (remove alias variables/equations as well as redundant equations,
  provide definite values to variables that have an infinite number of solutions if this makes sense,
  make state constraints structurally visible).
  Many equations of object-oriented models are linear Integer equations and can be pre-processed
  exactly to simplify the equations and to remove (consistently) redundant or
  overdetermined equations.
  
- Reduce the dimension of algebraic equation systems by tearing.

- Find equations that need to be differentiated one or more times (with the algorithm of Pantelides)
  in order that the DAE can be transformed to an ODE.
  
- Analytically differentiate the found equations.

- Statically select ODE states and transform to ODE form
  (hereby identifying linear equation systems that must be solved during simulation).
  

## Installation
 
The package is registered and is installed with (Julia >= 1.5 is required):

```julia
julia> ]add ModiaBase
```

It is recommended to also add the following packages, in order that all tests and examples can be executed in an own environment (`]test ModiaBase` works without adding these packages).

```julia
julia> ]add Unitful, Measurements, MonteCarloMeasurements, Distributions
```


## Main Developers

- Hilding Elmqvist, [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).

License: MIT (expat)
