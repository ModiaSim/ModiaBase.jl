# Transformation to ODE System

This section provides functions to **transform DAE to ODE systems**. In particular,

- to find equations that need to be differentiated, in order to transform a
  Differential Algebraic Equations (DAEs) agebraically to
  Ordinary Differential Equations (ODEs),

- to differentiate relevant equations analytically,

- to statically select ODE states and transform to ODE form.


## Main Functions

```@meta
CurrentModule = ModiaBase.BLTandPantelides
```

```@docs
pantelides!
```


```@meta
CurrentModule = ModiaBase.Differentiate
```

```@docs
derivative
```


```@meta
CurrentModule = ModiaBase
```

```@docs
StateSelectionFunctions
getSortedAndSolvedAST
```


## Utility Functions


```@meta
CurrentModule = ModiaBase
```

```@docs
StateCategory
ResidualCategory
LinearEquations
LinearEquationsIterator
EquationInfo
StateElementInfo
```
