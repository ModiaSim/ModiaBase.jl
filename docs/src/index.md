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

These functions are used by package [TinyModia](https://github.com/ModiaSim/TinyModia.jl),
but can also be utilized in another context.

Especially the following functionality is provided:

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

Transformation from a DAE to an ODE form is performed if no nonlinear-algebraic equations
appear and the ODE-states can be statically selected.

Array variables and array equations are kept (they are not "flattened" in single elements).
However, DAE forms that require to differentiate array equations, are not yet supported.

The following extensions are planned for the near future (a "few months"; most of the code is
available but must be adapted to ModiaBase):

- Full support of array equations.
- State- and time events.
- If transformation to an ODE is not possible with the algorithms above,
  transformation to a special index 1 DAE, that
  can be simulated with standard DAE solvers (such as Sundials IDA).


## Installation

The package is registered and is installed with (Julia >= 1.5 is required):

```julia
julia> ]add ModiaBase
```


It is recommended to also add the following packages, in order that all tests and examples can be executed:

```julia
julia> ]add Unitful, Measurements, MonteCarloMeasurements, Distributions
```

## Release Notes

### Version 0.7.1

- Due to version conflicts, added version 0.17 of DataStructures in compat.


### Version 0.7.0

- Initial version, based on code developed for Modia 0.6 and ModiaMath 0.6.


## Main developers

- Hilding Elmqvist (Mogram AB)
- Martin Otter ([DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
