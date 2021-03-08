module TestToStateSpace

using ModiaBase
using Test

@testset "\nTest toStateSpace.jl" begin

    #=
    Two inertias connected by an ideal gear
    
    der(phi1)  = w1
    der(phi2)  = w2
    J1*der(w1) = -tau1 + u
    J2*der(w2) = tau2
    0 = -phi1 + r*phi2
    0 = -r*tau1 + tau2
    
    Correct solution:
        w1 = der(phi1)
            tau2 = r*tau1
            J2*der(w2) = tau2 = r*tau1
            tau1 = J2/r*der(w2)
                = J2/r*der(w1)/r
                = J2/r^2*der(w1)
            J1*der(w1) = -tau1 + u
                    = -J2/r^2*der(w1) + u
            (J1 + J2/r^2)*der(w1) = u
        der(w1) = 1/(J1 + J2/r^2) * u
                = 0.3121019108280255*u
        y1 = tau1
        = J2/r^2*der(w1)
        = J2/r^2*0.3121019108280255*u
        = 0.06369426751592358*u
        y2 = phi1       
    
    x_name = [:phi1, :phi2, :w1, :w2, :tau1, :tau2]
    =#
    
    J1 =  3.0;
    J2 = 10.0;
    r  =  7.0;   # gear ratio
    
    E1 = [1 0  0  0 0 0
          0 1  0  0 0 0
          0 0 J1  0 0 0
          0 0  0 J2 0 0
          0 0  0  0 0 0
          0 0  0  0 0 0]
    
    A1 = [0 0 1 0  0 0
          0 0 0 1  0 0
          0 0 0 0 -1 0
          0 0 0 0  0 1
          -1 r 0 0  0 0
          0 0 0 0 -r 1]
    
    B1 = reshape([0, 0, 1.0, 0, 0, 0],6,1)
    
    C1 = [0   0 0 0 1.0 0;   # y1 = tau1
          1.0 0 0 0 0   0]   # y2 = phi1
    
    D1 = zeros(2,1)
    
    
    # First test
    (A2,B2,C2,D2,p2) = toStateSpace(E1, A1, B1, C1, D1)
    
    abstol = 1e-10
    k = 1/(J1 + J2/r^2)
    @test isapprox(A2, [0.0 0.0;1.0 0.0]           , atol=abstol)
    @test isapprox(B2, reshape([k,0],2,1)          , atol=abstol)
    @test isapprox(C2, [0.0 0.0; 0.0 1.0]          , atol=abstol)
    @test isapprox(D2, reshape([J2/r^2*k; 0.0],2,1), atol=abstol) 
    @test p2 == [3,1]
    
    # Second test
    (A3,B3,C3,D3,p3) = toStateSpace([1.0 0.0;0.0 1.0], A2, B2, C2, D2)
    @test isapprox(A3, A2, atol=abstol)
    @test isapprox(B3, B2, atol=abstol)
    @test isapprox(C3, C2, atol=abstol)
    @test isapprox(D3, D2, atol=abstol)    
    @test p3 == [1,2]
    
    # Third test
        #=
        # Find arbitrary regular matrix
           S = rand(6,6)
           while cond(S) > 1e6
              S = rand(6,6)
           end
          @show S
        =#
    S = [0.04657712210063569 0.5977856114171689 0.4786954916224406 0.04284939567378343 0.0695733063167907 0.43956659685818256; 
         0.5787577182369117 0.8137678461414808 0.4806941411619978 0.6472825897209098 0.5896925199894048 0.7869391012985769; 
         0.689055872370582 0.3725597060839918 0.03468538463077531 0.9068630223342913 0.6083998830118467 0.8483843026618747;
         0.22390058951722946 0.6851371355207307 0.5880543138208856 0.11322673490297497 0.33179865225658967 0.5333918928054828; 
         0.42423565831828003 0.6321345004194057 0.8419777975896683 0.3094611960849114 0.5643223914924458 0.4373555938996141; 
         0.4172735713046267 0.794126703767219 0.7967916197035996 0.17448224138354762 0.3026262810569218 0.9738148273956648]
         
    (A4,B4,C4,D4,p4) = toStateSpace(S*E1, S*A1, S*B1, C1, D1)
    @test isapprox(A4, [0.0 0.0;1.0 0.0]           , atol=abstol)
    @test isapprox(B4, reshape([k,0],2,1)          , atol=abstol)
    @test isapprox(C4, [0.0 0.0; 0.0 1.0]          , atol=abstol)
    @test isapprox(D4, reshape([J2/r^2*k; 0.0],2,1), atol=abstol) 
    @test p2 == [3,1]
end

end
