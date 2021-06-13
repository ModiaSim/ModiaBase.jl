"""
    module TestLinearIntegerEquations - Test ModiaBase/src/LinearIntegerEquations.jl
"""
module TestLinearIntegerEquations

using Test
using ModiaBase


println("... Test LinearIntegerEquations.jl")


@testset "Test LinearIntegerEquations.jl" begin

    @testset "... Test Voltage source and resistor without ground" begin
        println("\n    --- Test Voltage source and resistor without ground")
        #=            
            # Resistor
            1: R.v = R.p.v - R.n.v
            2: 0 = R.p.i + R.n.i
            3: R.v = R.R*R.p.i
            
            # Voltage source
            4: V.v = V.p.v - V.n.v
            5: 0   = V.p.i + V.n.i
            6: V.v = 10
            
            # Connect equations
            connect(V.p, R.p)
            connect(R.n, V.n)
            ->  7: V.p.v = R.p.v
                8: 0 = V.p.i + R.p.i
                9: R.n.v = V.n.v
            10: 0 = R.n.i + V.n.i
            
            # Variables:
            1: R.v
            2: R.p.v
            3: R.p.i
            4: R.n.v
            5: R.n.i
            6: V.v
            7: V.p.v
            8: V.p.i
            9: V.n.v
            10: V.n.i
        =#
        G = Vector{Int}[[1, 2, 4],
                        [3, 5],
                        [1, 3],
                        [6, 7, 9],
                        [8, 10],
                        [6],
                        [7, 2],
                        [8, 3],
                        [4, 9],
                        [5, 10]]
            
        eInt = Int[1,2,4,5,7,8,9,10]
        
        GcInt = Vector{Int}[[-1, 1, -1],
                            [1,1],
                            [-1,1,-1],
                            [1,1],
                            [-1,1],
                            [1,1],
                            [-1,1],
                            [1,1]]
                            
        Avar = fill(0,10)
        
        vNames = ["R.v",
                "R.p.v",
                "R.p.i",
                "R.n.v",
                "R.n.i",
                "V.v",
                "V.p.v",
                "V.p.i",
                "V.n.v",
                "V.n.i"]
                
        name(v) = vNames[v]
    
        (vEliminated, vProperty, nvArbitrary, redundantEquations) = simplifyLinearIntegerEquations!(G, eInt, GcInt, Avar)
        
        printSimplifiedLinearIntegerEquations(G, eInt, GcInt, vEliminated, vProperty, nvArbitrary, 
                                              redundantEquations, name, printTest=false)

        @test nvArbitrary == 1
        @test sort(vEliminated) == sort([9, 6, 2, 4, 10, 7, 8, 5])
        @test vProperty[[9, 6, 2, 4, 10, 7, 8, 5]] == [0, 1, 1, 0, 3, 1, -3, -3]
        @test redundantEquations == [10]
        @test sort(eInt)  == sort([2, 5, 7, 8, 9, 1, 4, 10])
        @test G[eInt] == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]
        @test GcInt == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]
      
    end
    
    
    
    @testset "... Test Voltage source and resistor with ground" begin
        println("\n    --- Test Voltage source and resistor with ground")
        #=  
            # Resistor
            1: R.v = R.p.v - R.n.v
            2: 0 = R.p.i + R.n.i
            3: R.v = R.R*R.p.i
            
            # Voltage source
            4: V.v = V.p.v - V.n.v
            5: 0   = V.p.i + V.n.i
            6: V.v = 10
            
            # Ground
            7: ground.p.v = 0
            
            # Connect equations
            connect(V.p, R.p)
            connect(R.n, V.n)
            connect(R.n, ground.p)
            ->   8: V.p.v = R.p.v
                 9: 0 = V.p.i + R.p.i
                10: R.n.v = V.n.v
                11: R.n.v = ground.p.v
                12: 0 = R.n.i + V.n.i + ground.p.i
            
            # Variables:
            1: R.v
            2: R.p.v
            3: R.p.i
            4: R.n.v
            5: R.n.i
            6: V.v
            7: V.p.v
            8: V.p.i
            9: V.n.v
            10: V.n.i
            11: ground.p.v
            12: ground.p.i
        =#
        G = Vector{Int}[[1, 2, 4],
                        [3, 5],
                        [1, 3],
                        [6, 7, 9],
                        [8, 10],
                        [6],
                        [11],
                        [7, 2],
                        [8, 3],
                        [4, 9],
                        [4,11],
                        [5, 10, 12]]
            
        eInt = Int[1,2,4,5,7,8,9,10,11,12]
        
        GcInt = Vector{Int}[[-1, 1, -1],
                            [1,1],
                            [-1,1,-1],
                            [1,1],
                            [1],
                            [-1,1],
                            [1,1],
                            [-1,1],
                            [-1,1],
                            [1,1,1]]
                            
        Avar = fill(0,12)
        
        vNames = ["R.v",
                "R.p.v",
                "R.p.i",
                "R.n.v",
                "R.n.i",
                "V.v",
                "V.p.v",
                "V.p.i",
                "V.n.v",
                "V.n.i",
                "ground.p.v",
                "ground.p.i"]
                
        name(v) = vNames[v]
    
        (vEliminated, vProperty, nvArbitrary, redundantEquations) = simplifyLinearIntegerEquations!(G, eInt, GcInt, Avar)
        
        printSimplifiedLinearIntegerEquations(G, eInt, GcInt, vEliminated, vProperty, nvArbitrary, 
                                              redundantEquations, name, printTest=false)

        @test nvArbitrary == 0
        @test sort(vEliminated) == sort([6, 12, 5, 10, 7, 2, 8, 9, 4, 11])
        @test vProperty[[6, 12, 5, 10, 7, 2, 8, 9, 4, 11]] == [1, 0, -3, 3, 1, 1, -3, 0, 0, 0]
        @test redundantEquations == Int64[]
        @test sort(eInt)  == sort([7, 11, 10, 5, 1, 8, 9, 2, 12, 4])
        @test G[eInt] == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]
        @test GcInt == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]
      
    end


    @testset "... Test Voltage source and capacitor without ground" begin
        println("\n    --- Test Voltage source and capacitor without ground")
        #=            
            # Capacitor
            1: C.v = C.p.v - C.n.v
            2: 0 = C.p.i + C.n.i
            3: C.C*der(C.v) = C.p.i
            
            # Voltage source
            4: V.v = V.p.v - V.n.v
            5: 0   = V.p.i + V.n.i
            6: V.v = 10
            
            # Connect equations
            connect(V.p, C.p)
            connect(C.n, V.n)
            ->  7: V.p.v = C.p.v
                8: 0 = V.p.i + C.p.i
                9: C.n.v = V.n.v
               10: 0 = C.n.i + V.n.i
            
            # Variables:
             1: C.v
             2: C.p.v
             3: C.p.i
             4: C.n.v
             5: C.n.i
             6: V.v
             7: V.p.v
             8: V.p.i
             9: V.n.v
            10: V.n.i
            11: der(C.v)
        =#
        G = Vector{Int}[[1, 2, 4],
                        [3, 5],
                        [11, 3],
                        [6, 7, 9],
                        [8, 10],
                        [6],
                        [7, 2],
                        [8, 3],
                        [4, 9],
                        [5, 10]]
            
        eInt = Int[1,2,4,5,7,8,9,10]
        
        GcInt = Vector{Int}[[-1, 1, -1],
                            [1,1],
                            [-1,1,-1],
                            [1,1],
                            [-1,1],
                            [1,1],
                            [-1,1],
                            [1,1]]
                    
        Avar = fill(0,10)
        pushfirst!(Avar,11)
        
        vNames = ["C.v",
                  "C.p.v",
                  "C.p.i",
                  "C.n.v",
                  "C.n.i",
                  "V.v",
                  "V.p.v",
                  "V.p.i",
                  "V.n.v",
                  "V.n.i",
                  "der(C.v)"]
                
        name(v) = vNames[v]
    
        (vEliminated, vProperty, nvArbitrary, redundantEquations) = simplifyLinearIntegerEquations!(G, eInt, GcInt, Avar)
        
        printSimplifiedLinearIntegerEquations(G, eInt, GcInt, vEliminated, vProperty, nvArbitrary, 
                                              redundantEquations, name, printTest=false)
        
        @test nvArbitrary == 1
        @test sort(vEliminated) == sort([9, 6, 2, 4, 10, 7, 8, 5])
        @test vProperty[[9, 6, 2, 4, 10, 7, 8, 5]] == [0, 1, 1, 0, 3, 1, -3, -3]
        @test redundantEquations == [10]
        @test sort(eInt)  == sort([2, 5, 7, 8, 9, 1, 4, 10])
        @test G[eInt] == Any[Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]
        @test GcInt == Any[Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]]
    end
    
    
    
    @testset "... Test two inductances in series with parallel resistors without ground" begin
        println("\n    --- Test two inductances in series with parallel resistors without ground")
        #=            
            # Inductor 1
            1: L1.v = L1.p.v - L1.n.v
            2: 0 = L1.p.i + L1.n.i
            3: L1.L*der(L1.p.i) = L1.v

            # Inductor 2
            4: L2.v = L2.p.v - L2.n.v
            5: 0 = L2.p.i + L2.n.i
            6: L2.L*der(L2.p.i) = L2.v
 
            # Resistor 1
            7: R1.v = R1.p.v - R1.n.v
            8:    0 = R1.p.i + R1.n.i
            9: R1.v = R1.R*R1.p.i
 
            # Resistor 2
            10: R2.v = R2.p.v - R2.n.v
            11:    0 = R2.p.i + R2.n.i
            12: R2.v = R2.R*R2.p.i
            
            # Voltage source
            13: V.v = V.p.v - V.n.v
            14: 0   = V.p.i + V.n.i
            15: V.v = 10
            
            # Connect equations
            connect(V.p, L1.p)
            connect(L1.n, R1.p)
            connect(L1.n, R2.p)
            connect(R1.n, L2.p)
            connect(R2.n, L2.p)
            connect(L2.n, V.n)
            ->  16: V.p.v = L1.p.v
                17: 0 = V.p.i + L1.p.i
                18: L1.n.v = R1.p.v
                19: L1.n.v = R2.p.v                
                20: 0 = L1.n.i + R1.p.i + R2.p.i
                21: R1.n.v = L2.p.v
                22: R2.n.v = L2.p.v                
                23: 0 = R1.n.i + R2.n.i + L2.p.i
                24: L2.n.v = V.n.v
                25: 0 = L2.n.i + V.n.i
                
            # Variables:
             1: L1.v
             2: L1.p.v
             3: L1.p.i
             4: L1.n.v
             5: L1.n.i
             6: L2.v
             7: L2.p.v
             8: L2.p.i
             9: L2.n.v
            10: L2.n.i
            11: R1.v
            12: R1.p.v
            13: R1.p.i
            14: R1.n.v
            15: R1.n.i
            16: R2.v
            17: R2.p.v
            18: R2.p.i
            19: R2.n.v
            20: R2.n.i               
            21: V.v
            22: V.p.v
            23: V.p.i
            24: V.n.v
            25: V.n.i
            26: der(L1.p.i)
            27: der(L2.p.i)     

    Solution of first version of LinearIntegerEquations.jl
            +++ Remove singularities
    Linear Integer equations:
       1: 0 = -L1.v + L1.p.v - L1.n.v
       2: 0 = L1.p.i + L1.n.i
       4: 0 = -L2.v + L2.p.v - L2.n.v
       5: 0 = L2.p.i + L2.n.i
       7: 0 = -R1.v + R1.p.v - R1.n.v
       8: 0 = R1.p.i + R1.n.i
      10: 0 = -R2.v + R2.p.v - R2.n.v
      11: 0 = R2.p.i + R2.n.i
      13: 0 = -V.v + V.p.v - V.n.v
      14: 0 = V.p.i + V.n.i
      16: 0 = -V.p.v + L1.p.v
      17: 0 = V.p.i + L1.p.i
      18: 0 = -L1.n.v + R1.p.v
      19: 0 = -L1.n.v + R2.p.v
      20: 0 = L1.n.i + R1.p.i + R2.p.i
      21: 0 = -R1.n.v + L2.p.v
      22: 0 = -R2.n.v + L2.p.v
      23: 0 = R1.n.i + R2.n.i + L2.p.i
      24: 0 = -L2.n.v + V.n.v
      25: 0 = L2.n.i + V.n.i
    Unknown variables:
       1: L1.v
       2: L1.p.v    (to be solved by equations)
       4: L1.n.v    (to be solved by equations)
       3: L1.p.i    (potential state)
       5: L1.n.i    (to be solved by equations)
       6: L2.v
       7: L2.p.v    (to be solved by equations)
       9: L2.n.v    (to be solved by equations)
       8: L2.p.i    (potential state)
      10: L2.n.i    (to be solved by equations)
      11: R1.v
      12: R1.p.v    (to be solved by equations)
      14: R1.n.v    (to be solved by equations)
      13: R1.p.i
      15: R1.n.i    (to be solved by equations)
      16: R2.v
      17: R2.p.v    (to be solved by equations)
      19: R2.n.v    (to be solved by equations)
      18: R2.p.i
      20: R2.n.i    (to be solved by equations)
      21: V.v
      22: V.p.v    (to be solved by equations)
      24: V.n.v    (to be solved by equations)
      23: V.p.i    (to be solved by equations)
      25: V.n.i    (to be solved by equations)

    After first transformation to trapezoidal form (eliminate variables that must be solved):
       2: 0 = L1.n.i + L1.p.i
       5: 0 = L2.n.i + L2.p.i
       8: 0 = R1.n.i + R1.p.i
      11: 0 = R2.n.i + R2.p.i
      14: 0 = V.p.i + V.n.i
      16: 0 = -V.p.v + L1.p.v
      17: 0 = -V.n.i + L1.p.i
      18: 0 = -L1.n.v + R1.p.v
      19: 0 = R2.p.v - R1.p.v
      21: 0 = -R1.n.v + L2.p.v
      22: 0 = -R2.n.v + L2.p.v
      24: 0 = -L2.n.v + V.n.v
       1: 0 = L1.p.v - L1.v - R1.p.v
      13: 0 = -V.n.v - V.v + L1.v + R1.p.v
       7: 0 = -R1.p.v + R1.v + L2.p.v
    ---------- rk1
      20: 0 = R1.p.i + R2.p.i - L1.p.i
       4: 0 = L2.v - V.v + L1.v + R1.v
      23: 0 = L2.p.i - R1.p.i - R2.p.i
      10: 0 = R2.v - R1.v
      25: 0 = -L2.p.i + L1.p.i

    After second transformation to trapezoidal form (ignore potential states):
       2: 0 = L1.n.i + L1.p.i
       5: 0 = L2.n.i + L2.p.i
       8: 0 = R1.n.i + R1.p.i
      11: 0 = R2.n.i + R2.p.i
      14: 0 = V.p.i + V.n.i
      16: 0 = -V.p.v + L1.p.v
      17: 0 = -V.n.i + L1.p.i
      18: 0 = -L1.n.v + R1.p.v
      19: 0 = R2.p.v - R1.p.v
      21: 0 = -R1.n.v + L2.p.v
      22: 0 = -R2.n.v + L2.p.v
      24: 0 = -L2.n.v + V.n.v
       1: 0 = L1.p.v - L1.v - R1.p.v
      13: 0 = -V.n.v - V.v + L1.v + R1.p.v
       7: 0 = -R1.p.v + R1.v + L2.p.v
    ---------- rk1
      10: 0 = R2.v - R1.v
       4: 0 = L2.v - V.v + L1.v + R1.v
      23: 0 = -R1.p.i + L2.p.i - R2.p.i
    ---------- rk2
      20: 0 = L1.p.i - L2.p.i
      25: 0 = -L2.p.i + L1.p.i

    After third transformation to trapezoidal form (eliminate potential states):
       2: 0 = L1.n.i + L1.p.i
       5: 0 = L2.n.i + L2.p.i
       8: 0 = R1.n.i + R1.p.i
      11: 0 = R2.n.i + R2.p.i
      14: 0 = V.p.i + V.n.i
      16: 0 = -V.p.v + L1.p.v
      17: 0 = -V.n.i + L1.p.i
      18: 0 = -L1.n.v + R1.p.v
      19: 0 = R2.p.v - R1.p.v
      21: 0 = -R1.n.v + L2.p.v
      22: 0 = -R2.n.v + L2.p.v
      24: 0 = -L2.n.v + V.n.v
       1: 0 = L1.p.v - L1.v - R1.p.v
      13: 0 = -V.n.v - V.v + L1.v + R1.p.v
       7: 0 = -R1.p.v + R1.v + L2.p.v
    ---------- rk1
      10: 0 = R2.v - R1.v
       4: 0 = L2.v - V.v + L1.v + R1.v
      23: 0 = -R1.p.i + L2.p.i - R2.p.i
    ---------- rk2
      20: 0 = L1.p.i - L2.p.i
    ---------- rk3

    After alias elimination:
       1: 0 = L1.p.v - L1.v - R1.v
      13: 0 = -V.n.v - V.v + L1.v + R1.v
    ---------- rk1
       4: 0 = L2.v - V.v + L1.v + R1.v
      23: 0 = -R1.p.i + L2.p.i - R2.p.i
    ---------- rk2
      20: 0 = L1.p.i - L2.p.i
    ---------- rk3

    Final, simplified equations:
       1: 0 = L1.p.v - L1.v - R1.v
      13: 0 = -V.v - L1.p.v + V.n.v
    ---------- rk1
      23: 0 = -R1.p.i + L2.p.i - R2.p.i
    ---------- rk2
      20: 0 = L1.p.i - L2.p.i
    ---------- rk3

    Variables that can be arbitrarily set and have been set to zero:
       7: L2.p.v = 0

    Variables that have been eliminated:
      16: R2.v = R1.v
      12: R1.p.v = R1.v
       9: L2.n.v = V.n.v
      19: R2.n.v = 0
      14: R1.n.v = 0
      17: R2.p.v = R1.v
       4: L1.n.v = R1.v
      25: V.n.i = L1.p.i
      22: V.p.v = L1.p.v
      23: V.p.i = -L1.p.i
      20: R2.n.i = -R2.p.i
      15: R1.n.i = -R1.p.i
      10: L2.n.i = -L2.p.i
       5: L1.n.i = -L1.p.i
       6: L2.v = -V.n.v

    Redundant equations that have been removed:
      25

    Remaining transformed linear Integer equations:
       1: 0 = L1.p.v - L1.v - R1.v
      13: 0 = -V.v + L1.p.v - V.n.v
      23: 0 = -R1.p.i + L2.p.i - R2.p.i
      20: 0 = L1.p.i - L2.p.i
        =#
        G = Vector{Int}[[1, 2, 4],
                        [3, 5],
                        [26, 1],
                        
                        [6, 7, 9],
                        [8,10],
                        [27, 6],
                        
                        [11, 12, 14],
                        [13, 15],
                        [11, 13],
                        
                        [16, 17, 19],
                        [18, 20],
                        [16, 18],
                        
                        [21, 22, 24],
                        [23, 25],
                        [21],
                        
                        [22, 2],
                        [23, 3],
                        
                        [4, 12],
                        [4, 17],
                        [5, 13, 18],
                        
                        [14, 7],
                        [19, 7],
                        [15, 20, 8],
                        
                        [9, 24],
                        [10, 25]]
            
        eInt = Int[1,2,4,5,7,8,10,11,13,14,16,17,18,19,20,21,22,23,24,25]
        
        GcInt = Vector{Int}[[-1, 1, -1],
                            [1,1],
                            [-1,1,-1],
                            [1,1],
                            [-1,1,-1],
                            [1,1],
                            [-1,1,-1],
                            [1,1],
                            [-1,1,-1],
                            [1,1],
                            [-1,1],
                            [1,1],
                            [-1,1],
                            [-1,1],
                            [1,1,1],
                            [-1,1],
                            [-1,1],
                            [1,1,1],
                            [-1,1],
                            [1,1]]
                    
        Avar = fill(0,27)
        Avar[3] = 26
        Avar[8] = 27
        
        vNames = ["L1.v", 
                  "L1.p.v", 
                  "L1.p.i",
                  "L1.n.v",
                  "L1.n.i",
                  "L2.v",
                  "L2.p.v",
                  "L2.p.i",
                  "L2.n.v",
                  "L2.n.i",
                  "R1.v",
                  "R1.p.v",
                  "R1.p.i",
                  "R1.n.v",
                  "R1.n.i",
                  "R2.v",
                  "R2.p.v",
                  "R2.p.i",
                  "R2.n.v",
                  "R2.n.i",               
                  "V.v",
                  "V.p.v",
                  "V.p.i",
                  "V.n.v",
                  "V.n.i",
                  "der(L1.p.i)",
                  "der(L2.p.i)"] 
    
        var_name(v) = vNames[v]
        (vEliminated, vProperty, nvArbitrary, redundantEquations) = simplifyLinearIntegerEquations!(G, eInt, GcInt, Avar; log=true, var_name = var_name)
        
        printSimplifiedLinearIntegerEquations(G, eInt, GcInt, vEliminated, vProperty, nvArbitrary, redundantEquations, var_name, printTest=false)

        #= previously:
            @test vEliminated == [7, 16, 12, 9, 19, 14, 17, 4, 25, 22, 23, 20, 15, 10, 5, 6]
            @test Property[vEliminated] = [0, 11, 11, 24, 0, 0, 11, 11, 3, 2, -3, -18, -13, -8, -3, -24]
            @test eInt == [2, 5, 8, 11, 14, 16, 17, 18, 19, 21, 22, 24, 1, 13, 7, 10, 4, 23, 20, 25]
            @test G[eInt] == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], [2, 1, 11], [21, 2, 24], Int64[], Int64[], Int64[], [13, 8, 18], [3, 8], Int64[]]
            @test GcInt == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], [1, -1, -1], [-1, 1, -1], Int64[], Int64[], Int64[], [-1, 1, -1], [1, -1], Int64[]]     
        =#
        
        @test nvArbitrary == 1
        @test vEliminated == [2, 25, 23, 20, 15, 10, 5, 16, 24, 9, 19, 14, 17, 4, 22, 1]
        @test vProperty[vEliminated] == [0, 3, -3, -18, -13, -8, -3, 11, -21, -21, 7, 7, 12, 12, 0, -12]
        @test redundantEquations == [25]
        @test eInt  == [16, 2, 18, 5, 19, 8, 21, 11, 22, 14, 24, 17, 4, 7, 20, 13, 10, 23, 1, 25]
        @test G[eInt] == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], [7, 6, 21], [11, 12, 7], [13, 18, 3], Int64[], Int64[], [8, 3], Int64[], Int64[]]
        @test GcInt == [Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], Int64[], [1, -1, 1], [-1, 1, -1], [1, 1, -1], Int64[], Int64[], [1, -1], Int64[], Int64[]]
    end
end

end
