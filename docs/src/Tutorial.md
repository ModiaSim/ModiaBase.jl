# Tutorial

This chapter contains a short tutorial about the data structures and
functions provided by package ModiaBase.


## 1. Regular DAEs (Index Zero DAEs)

In this subsection functions are demonstrated that can be used to
transform regular DAEs to ODEs.

The transformations are explained with the
following simple electrical circuit (a low pass filter where the
voltage source is modelled with an inner resistor):

![Low Pass Filter](../resources/images/LowPassFilter.png)



### 1.1 Bi-Partite Graph

After simple equations have been removed, this circuit can be
described with 6 equations. The structure of the equations is
displayed in the next figure:

![Incidence Matrix of Low Pass Filter](../resources/images/LowPassFilterReduced_IncidenceMatrix.png)

- Every **column** corresponds to one time-varying **variable**.
  Parameters, so variables with constant values, are not shown.

- Every **row** corresponds to one **equation**.

- A cell is marked (here in *blue*), if a time-varying variable is
  present in one equation. Variables that are appearing differentiated,
  such as `C.v`, are not marked because in a first analysis phase, these potential
  state variables are treated as known.

The matrix above is called the **incidence matrix** or the
**bi-partite graph** of the circuit. In ModiaBase, this matrix is represented
as vector `G` of integer vectors:

```julia
# Bi-partite graph of low pass filter
G = [ [1,2,4],  # equation 1 depends on variables 1,2,4
      [1,7],
      [3,4],
      [3,7],
      [6,7],
      [2] ]
```

This can be also made more explicit (and a bit more efficient storage
by defining the incidence matrix as):


```julia
# Bi-partite graph of low pass filter
G = Vector{Int}[ [1,2,4],
                 [1,7],
                 [3,4],
                 [3,7],
                 [6,7],
                 [2] ]
```


### 1.2 Assignment

In a first step, an assignment is made (also called matching), to associate
one variable uniquely with one equation:

![Matched IIncidence Matrix of Low Pass Filter](../resources/images/LowPassFilterReduced_Matching.png)

- *Red* marks show the assigned variables.
- *Blue* marks show if a variable is part of the respective equation

The assignment is computed with function [`ModiaBase.matching`](@ref)
returning a vector **assign**:

```julia
using ModiaBase
M          = 7  # Number of variables
vActive    = fill(true,M)
vActive[5] = false    # state C.v is known
assign     = matching(G, M, vActive)

# assign = [2,6,3,1,0,5,4]
```

The meaning of vector `assign` is that

- Variable 1 is solved from equation 2,
- Variable 2 is solved from equation 6,
- etc.


### 1.3 Block Lower Triangular transformation

In a second step, equations are sorted and algebraic loops determined:

![Incidence Matrix of sorted equations of Low Pass Filter](../resources/images/LowPassFilterReduced_BLT.png)

- *Red* marks show the assigned variables.
- *Blue* marks show if a variable is part of the respective equation
- A *grey* area marks an algebraic loop.

The sorting is computed with function [`ModiaBase.BLT`](@ref):

```julia
using ModiaBase
blt = BLT(G, assign)

#=
    blt = [ [6],
           [3,4,2,1],
           [5] ]
=#
```

The meaning is for example that the second BLT block consists of
equations 3,4,2,1 and these equations form an algebraic loop.



## 2. Singular DAEs (Higher Index DAEs)

xxx
