# Data Structures

In this chapter the **basic data structures** are summarized and
shortly described that are used in package ModiaBase.


## 1. Bi-partite Graph

The *bi-partite Graph* `G` of a DAE system defines the functional dependency of the equations ``e_i`` from
time-varying variables ``v_j``. `G` is also a sparse representation of the incidence matrix of the DAE system.
Example:

```julia
# Bi-partite graph of low pass filter
G = Vector{Int}[ [1,2,4], # equation 1 depends on variables 1,2,4
                 [1,7],   # equation 2 depends on variables 1,7
                 [3,4],
                 [3,7],
                 [6,7],
                 [2] ]
```

In ModiaBase, only potential states and unknown variables are described with `G`.
Dependency on parameters (= constant quantities) is not shown.

The number of variables is usually larger as the number of equations, because both
the potential states and their derivatives are stored in `G`.


## 2. Assignment Vector

The *assignment vector* `assign` defines for all equations and variables the unique variable ``v_j`` that
is solved from equation ``e_i``. Example:

```julia
assign = [2,6,3,1,0,5,4]  # Variable 1 is solved from equation 2
                          # Variable 2 is solved from equation 6
                          #  ...
                          # Variable 5 is not solved from any equation
```

The *inverted assignment vector* `invAssign` defines the unique equation ``e_i`` that
is solved for the unique variable ``v_j``. This vector can be directly computed
from vector `assign`. Example:


```julia
(invAssign, unAssigned) = invertAssign(assign)

# invAssign = [4,1,3,7,6,2,0]   # Equation 1 is solved for variable 4
                                # Equation 2 is solved for variable 1
                                # ...
                                # Equation 7 is not solved for any variable
```


## 3. Block Lower Triangular Form

The *Block Lower Triangular form* `blt` of an equation system describes the
sorted set of equations, in order to solve for the unknown variables.
With vector `invAssign` (see subsection 2 above) the information is provided
for which variable the respective equation is solved.

Example:

```julia
blt = [ [6],        # Solve first equation 6
        [3,4,2,1],  # Afterwards solve equations 3,4,2,1 (they form an algebraic loop)
        [5] ]       # Finally solve equation 5

invAssign = [4, 1, 3, 7, 6, 2, 0]
```

The meaning is:

1. Equation 6 is solved for variable 2.
2. Equations 3,4,2,1 are solved simultaneously for variables 3, 7, 1, 4.
3. Equation 5 is solved for variable 6.


## 4. Variable Association Vector

The derivative relationship between variables is described with
the *variable association vector* `Avar` and its
inverted vector `invAvar`:

```math
\begin{aligned}
    Avar_j &= \left\{ \begin{array}{rl}
                        k & \text{\textbf{if}}~ \dot{v}_j \equiv v_k \\
                        0 & \text{\textbf{if}}~ v_j \text{ is not a differentiated variable}
                      \end{array}\right. \\
    invAvar_j &= \left\{ \begin{array}{rl}
                            k & \text{\textbf{if}}~ v_j \equiv \dot{v}_k \\
                            0 & \text{\textbf{if}}~ v_j \text{ is not a differentiated variable}
                         \end{array}\right.
\end{aligned}
```

Example:

The following derivative relationships between variables `v1,v2,v3,v4,v5`

```julia
   1. v1
   2. v2 = der(v1)
   3. v3 = der(v2)
   4. v4
   5. v5 = der(v4)
```

are expressed by the following variable association vector and its inverted form:

```julia
   Avar    = [2,3,0,5,0]   # The derivative of variable 1 is variable 2
                           # The derivative of variable 2 is variable 3
                           #  ...

   invAvar = [0,1,2,0,4]   # Variable 1 is not a derivative
                           # Variable 2 is the derivative of variable 1
                           #   ...
                           # Variable 5 is the derivative of variable 4
```


## 5. Equation Association Vector

The derivative relationship between equations is described with
the *equation association vector* `Bequ` and its
inverted vector `invBequ`:

```math
\begin{aligned}
    Bequ_i &= \left\{ \begin{array}{rl}
                        k & \text{\textbf{if}}~ \dot{e}_i \equiv e_k \\
                        0 & \text{\textbf{if}}~ \dot{e}_i \text{ does not exist}
                      \end{array}\right. \\
    invBequ_i &= \left\{ \begin{array}{rl}
                            k & \text{\textbf{if}}~ e_i \equiv \dot{e}_k \\
                            0 & \text{\textbf{if}}~ _i \text{ is not a differentiated equation}
                         \end{array}\right.
\end{aligned}
```

Example:

The following derivative relationships between equations `e1,e2,e3,e4,e5`

```julia
   1. e1
   2. e2 = der(e1)
   3. e3 = der(e2)
   4. e4
   5. e5 = der(e4)
```

are expressed by the following equation association vector and its inverted form:

```julia
   Bvar    = [2,3,0,5,0]   # The derivative of equation 1 is equation 2
                           # The derivative of equation 2 is equation 3
                           #  ...

   invBvar = [0,1,2,0,4]   # Equation 1 is not a differentiated equation
                           # Equation 2 is the derivative of equation 1
                           #   ...
                           # Equation 5 is the derivative of equation 4
```
