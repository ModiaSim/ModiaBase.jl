# Equation Reduction

This section provides functions to **reduce the dimensions of
systems of equations**.



## Linear Integer Equations

Many equations of object-oriented models are **linear Integer equations**
(such as all equations deduced from connections, or, say defining a potential difference)
and can be pre-processed exactly to simplify the equations, for example elimination of
alias variables, or variables that are identically to zero). Furthermore,
(consistently) redundant or (consistently) overdetermined equations can be
removed. Finally, hidden state constraints can be made explicit in order
that a structural algorithm (such as the algorithm of Pantelides) can process state constraints.

```@meta
CurrentModule = ModiaBase
```

```@docs
simplifyLinearIntegerEquations!
printSimplifiedLinearIntegerEquations
```


## Algebraic Systems

**Algebraic equation systems** are reduced by selecting a subset of the
variables as iteration variables and computing the remaining variables in
a forward sequence.

```@meta
CurrentModule = ModiaBase
```

```@docs
TearingSetup
tearEquations!
```
