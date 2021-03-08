# License for this file: MIT (expat)
# Copyright 2021, DLR Institute of System Dynamics and Control
# Author: Martin Otter, DLR-SR


using LinearAlgebra
export toStateSpace

"""
    (Ar,Br,Cr,Dr,pr) = toStateSpace(E,A,B,C,D)
    
Transform the descriptor system

```julia
E*dx/dt = A*x + B*u
      y = C*x + D*u
```

where `E` may be singular, into state space form

```julia
dxr/dt = Ar*xr + Br*u
     y = Cr*xr + Dr*u
    xr = x[pr]
```

The elements of `xr` are a subset of the elements of `x`
characterized by the index vector `pr`.

`E, A, B, C, D` must be of type `Matrix{T<:AbstractFloat}` of
appropriate dimensions.


The function triggers an error in the following cases:

- The system has only algebraic equations (for example, if `E=0`).

- The system has no unique solution 
  (that is, it has none or an infinite number of solutions;
   for example, if it has redundant equations).

- A transformation to state space form with `xr` a subset of `x` would introduce a term
  `B2*der(u)`, that is the derivative of `u` would be required.
  
  
# Main developer

[Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/), 
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)  


# Algorithm

The algorithm is an extended version published in:

Martin Otter (2001): Kapitel 20.8 - Singuläre Deskriptorsysteme. Section in book:
   Dierk Schröder: Elektrische Antriebe - Regelung von Antriebssystemen, 2. Auflage.
"""
function toStateSpace(E::Matrix{T}, A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}) where {T<:AbstractFloat}
    #=
    1. QR decomposition of E
       E*der(x) = Q1*R*der(x)[p1] = A[:,p1]*x[p1] + B*u
       R*der(x)[p1] = transpose(Q1)*A[:,p1]*x[p1] + transpose(Q1)*B*u
       [R1
        0]*der(x)[p1] = A1*x[p1] + B1*u
                        A1 = transpose(Q1)*A[:,p1]
                        B1 = transpose(Q1)*B
                        n1 = size(R1,1)
                        x1 = x[p1]
       [R1
        0]*der(x1) = [A11
                      A21]*x1 + [B11
                                 B21]*u
                        A11 = A1[1:n1,:]
                        A21 = A1[n1+1:end,:]
                        B11 = B1[1:n1,:]
                        B21 = B1[n1+1:end,:]
                        
       R1*der(x1) = A11*x1 + B11*u
                0 = A21*x1 + B21*u
                          
    2. QR decomposition of A21
       0 = QQ*RR*x1[p2] + B21*u
           x2 = x1[p2]
              = x[p1][p2]
              = x[p1[p2]]
       0 = RR*x2 + transpose(QQ)*B21*u
       0 = [RR1 RR2]*[x12;x22] + B3*u
           n2 = size(RR1,2)
           B3 = transpose(QQ)*B21
           p3 = p1[p2]
           x12 = x2[1:n2]
               = x[p3[1:n2]]
           x22 = x2[n2+1:end]
               = x[p3[n2+1:end]]
           RR1 = RR[:,1:n2]
           RR2 = RR[:,n2+1:end]
           
       0 = RR1*x12 + RR2*x22 + B3*u    # RR1 must be non-singular
       x12 = RR1\(-RR2*x22 - B3*u)
           = A4*x22 + B4*u
             A4 = RR1\(-RR2)
             B4 = RR1\(-B3)
       
    3. Elimination of x12
       R1*der(x1) = A11*x1 + B11*u
       R1[:,p2]*der(x1)[p2] = A11[:,p2]*x1[p2] + B11*u
       R1[:,p2]*der(x2) = A11[:,p2]*x2 + B11*u
       R21*der(x12) + R22*der(x22) = A31*x12 + A32*x22 + B11*u
               R21 = R1[:,p2[1:n2]]
               R22 = R1[:,p2[n2+1:end]]
               A31 = A11[:,p2[1:n2]]
               A32 = A11[:,p2[n2+1:end]]

       R21*(A4*der(x22) + B4*der(u)) + R22*der(x22) = A31*(A4*x22 + B4*u) + A32*x22 + B11*u 
       (R21*A4 + R22)*der(x22) = (A31*A4 + A32)*x22 + (A31*B4 + B11)*u
                 Requirement: R21*B4 = 0
       Enew*der(xnew) = Anew*xnew + Bnew*u
           Enew = R21*A4 + R22
           Anew = A31*A4 + A32
           Bnew = A31*B4 + B11
           xnew = x22 = x[p3[n2+1:end]]
                      = x[pnew]
           pnew = p3[n2+1:end]

       y = C*x + D*u
         = C[:,p3]*x[p3] + D*u
         = [C11 C12]*[x12;x22] + D*u
           C11 = C[:,p3][:,1:n2]
           C12 = C[:,p3][:,n2+1:end]
       y = C11*x12 + C12*x22 + D*u
         = C11*(A4*x22 + B4*u) + C12*x22 + D*u
         = (C11*A4 + C12)*x22 + (C11*B4 + D)*u
         
       y = Cnew*xnew + Dnew*u
           Cnew = (C11*A4 + C12)
           Dnew = (C11*B4 + D)
    =#
    @assert(size(E,1) == size(E,2))
    @assert(size(A,1) == size(E,1))
    @assert(size(A,2) == size(A,1))
    @assert(size(B,1) == size(A,1))
    @assert(size(C,2) == size(A,1))
    @assert(size(D,1) == size(C,1))
    @assert(size(D,2) == size(B,2))
    
    E = deepcopy(E)
    A = deepcopy(A)
    B = deepcopy(B)
    C = deepcopy(C)
    D = deepcopy(D)
    p = 1:size(A,1)
    
    while true
        # 1. QR decomposition of E
        nx = size(E,1)       
        F  = qr(E, Val(true))
        Q1 = F.Q
        R  = F.R
        p1 = F.p 
        
        # Determine the lower part of R that is zero
        n1  = 0
        tol = eps(T)*nx*max(norm(A,Inf), norm(E,Inf))
        for i = nx:-1:1
            if abs(R[i,i]) > tol
                n1 = i
                break
            end
        end

        if n1 < 1
            # Temporarily trigger an error. However, would be possible 
            # to completely eliminate x and then compute Cr, Dr.
            error("toStateSpace(E,A,B,C,D): E*der(x) = A*x + B*u; y = C*x + D*u cannot be transformed\n",
                  "to state space form, because the system has only algebraic equations!")
        end
                  
        # Transform A and B
        A1 = transpose(Q1)*A[:,p1]
        B1 = transpose(Q1)*B
    
        # If R is regular, transform directly to state space form
        if n1 == nx
            # R*der(x)[p1] = A1*x[p1] + B1*u
            #            y = C*x + D*u
            #              = C[:,p1]*x[p1] + D*u
            Ar = R\A1
            Br = R\B1
            Cr = C[:,p1]
            pr = p[p1]
            return (Ar,Br,Cr,D,pr)
        end
        
        # Current status
        #   R1*der(x1) = A11*x1 + B11*u
        #            0 = A21*x1 + B21*u
        #                A11 = A1[1:n1,:]
        #                A21 = A1[n1+1:end,:]
        #                B11 = B1[1:n1,:]
        #                B21 = B1[n1+1:end,:]
        #
        # QR decomposition of A21
        R1  = R[1:n1,:]
        A21 = A1[n1+1:nx,:]
        B11 = B1[1:n1,:]
        B21 = B1[n1+1:nx,:]
        n2  = nx - n1
        F   = qr(A21, Val(true))
        QQ  = F.Q
        RR  = F.R
        p2  = F.p
        if abs(RR[n2,n2]) <= tol
            error("toStateSpace(E,A,B): E*der(x) = A*x + B*u; y = C*x + D*u cannot be transformed\n",
                  "to state space form, because the system has no unique solution!")        
        end
        
        # Current status
        #   R21*der(x12) + R22*der(x22) = A31*x12 + A32*x22 + B11*u
        #                             0 = RR1*x12 + RR2*x22 + B3*u 
        #       R21 = R1[:,p2[1:n2]]
        #       R22 = R1[:,p2[n2+1:end]]
        #       A31 = A11[:,p2[1:n2]]
        #       A32 = A11[:,p2[n2+1:end]]
        #       RR1 = RR[:,1:n2]
        #       RR2 = RR[:,n2+1:end]
        p12 = p2[1:n2]
        p22 = p2[n2+1:nx]
        R21 = R1[:,p12]
        R22 = R1[:,p22]
        A31 = A1[1:n1,p12]
        A32 = A1[1:n1,p22]
        B3  = transpose(QQ)*B21        
        RR1 = RR[:,1:n2]
        RR2 = RR[:,n2+1:nx]
        
        # Solve for x12
        #   x12 = RR1\(-RR2*x22 - B3*u)
        #       = A4*x22 + B4*u
        A4 = RR1\(-RR2)
        B4 = RR1\(-B3)          
        
        # Check that transformation to state space form is possible
        if norm(R21*B4, Inf) > tol
            error("toStateSpace(E,A,B): E*der(x) = A*x + B*u; y = C*x + D*u cannot be transformed\n",
                  "to state space form der(xr) = Ar*xr + Br*u; y = Cr*xr + Dr*u with xr a subset of x,\n",
                  "because the derivative of u is needed for the transformed system.")   
        end

        # Elimination of x12
        p3  = p1[p2] 
        C11 = C[:,p3[1:n2]]
        C12 = C[:,p3[n2+1:nx]]        
        E   = R21*A4 + R22
        A   = A31*A4 + A32
        B   = A31*B4 + B11
        C   = C11*A4 + C12
        D   = C11*B4 + D
        p   = p[ p3[n2+1:nx] ]          
    end
end
