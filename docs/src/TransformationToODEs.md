# Transformation to ODE System

This section provides utility functions to **transform DAE to ODE systems**. In particular,

- to find equations that need to be differentiated, in order to transform a
  Differential Algebraic Equations (DAEs) agebraically to
  Ordinary Differential Equations (ODEs),

- to differentiate relevant equations analytically,


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

