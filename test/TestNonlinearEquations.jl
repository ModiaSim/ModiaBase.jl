"""
    module TestNonlinearEquations

Module to test ModiaBase/src/NonlinearEquations.jl

"""
module TestNonlinearEquations

using Test
using ModiaBase.NonlinearEquations
using LinearAlgebra


println("... Test NonlinearEquations.jl")

@testset "\nTest NonlinearEquations.jl" begin

    loglevel = info

    # Linear underdetermined test
    @testset "... Test linear problem:" begin

        A = [1 2 3; 4 5 6]
        xref = [-1; -2; 5]
        b = A*xref
        m = size(A, 1)
        n = size(A, 2)

        # Test function
        function fx!(fx::Vector, x::Vector)
            @assert(length(x) == n)
            @assert(length(fx) == m)
            res = A*x - b
            for i in 1:m
                fx[i] = res[i]
            end
            return true
        end

        # Solver parameters
        tol = 1e-6
        maxiter = 50
        nonlinearity = mild
        restricted = false
        quasi = true
        forwardjac = false

        if loglevel >= info
            println("Linear problem:")
        end
        x = zeros(Float64, n)
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        @test isapprox(x, [-2.33333333333333, 0.666666666666666, 3.66666666666666])

    end


    # Determined test from https://www.zib.de/codelib/NewtonLib/
    @testset "... Test Chebyquad problem of dimensions 2..7" begin

        # Test function
        function fx!(fx::Vector, x::Vector)
            n = length(x)
            for i in range(1, n-1, step=2)
                fx[i] = 0.0
                fx[i+1] = n / ((i + 1.0)^2 - 1.0)
            end
            if isodd(n)
                fx[n] = 0.0
            end
            for l = 1:n
                factt = 4.0 * x[l] - 2.0
                ti2 = 1.0
                ti1 = 0.5 * factt
                fx[1] = ti1 + fx[1]
                for i = 2:n
                    ti = factt * ti1 - ti2
                    fx[i] = ti + fx[i]
                    ti2 = ti1
                    ti1 = ti
                end
            end
            return true
        end

        # Test parameters
        nn = 7  # max. Chebyquad problem dimension

        # Checks
        @assert(nn >= 2)

        # Solver parameters
        tol = 1e-6
        maxiter = 50
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = false

        for n = 2:nn

            if n == 6
                continue
            end

            if loglevel >= info
                println("Chebyquad problem n = $n:")
            end
            x = zeros(Float64, n)
            y = zeros(Float64, n)
            fscale = zeros(Float64, n)
            for i = 1:n
                x[i] = i / (n + 1)
                fscale[i] = 1.0
            end
            if loglevel == debug
                println("* start vector = $x")
                println("* scale vector = $fscale")
            end

            converged = solveNonlinearEquations!(fx!, n, x, fscale; nonlinearity=high, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

            if loglevel == debug
                println("* solution = $x")
                fx!(y, x)
                println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
            end

            if n == 2
                @test isapprox(x, [0.2113248654051863, 0.7886751345948136])
            elseif n == 3
                @test isapprox(x, [0.14644660948566998, 0.4999999999999999, 0.8535533905143302])
            elseif n == 4
                @test isapprox(x, [0.10267275999827528, 0.4062037728007276, 0.5937962271992724, 0.8973272400017247])
            elseif n == 5
                @test isapprox(x, [0.08375125622814243, 0.31272929406063793, 0.5, 0.6872707059393621, 0.9162487437718575])
            elseif n == 6
                @test converged == false
            elseif n == 7
                @test isapprox(x, [0.05806914954476061, 0.23517161265728312, 0.33804409443504935, 0.5, 0.6619559055649507, 0.764828387342717, 0.9419308504552394])
            elseif n == 8
                @test converged == false
            elseif n == 9
                @test converged == false
            end

        end
    end


    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/2by2.jl
    # and https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/already_converged.jl
    @testset "... Test 2by2 problem" begin

        # Test function
        m = 2
        n = 2
        function f_2by2!(F, x)
            F[1] = (x[1]+3)*(x[2]^3-7)+18
            F[2] = sin(x[2]*exp(x[1])-1)
            return true
        end

        # Solver parameters
        tol = 1e-7
        maxiter = 50
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = false

        for ii in 1:2
            if loglevel >= info
                println("2by2 problem $ii:")
            end
            if ii == 1
                x = [-0.5, 1.4]
            else
                x = [0.0, 1.0]  # already converged
            end
            y = zeros(Float64, m)
            fscale = ones(Float64, n)
            if loglevel == debug
                println("* start vector = $x")
                println("* scale vector = $fscale")
            end

            converged = solveNonlinearEquations!(f_2by2!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

            if loglevel == debug
                println("* solution = $x")
                f_2by2!(y, x)
                println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
            end

            @test isapprox(x, [0.0, 1.0], atol=100*tol)
        end

    end


    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
    @testset "... Test MINPACK Powell singular problem" begin

        # Test function
        m = 4
        n = 4
        function fx!(fvec, x)
            fvec[1] = x[1] + 10x[2]
            fvec[2] = sqrt(5)*(x[3] - x[4])
            fvec[3] = (x[2] - 2x[3])^2
            fvec[4] = sqrt(10)*(x[1] - x[4])^2
            return true
        end

        # Initial guess
        x = [3.0, -1.0, 0.0, 1.0]

        # Solver parameters
        tol = 1e-8
        maxiter = 50
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = false

        if loglevel >= info
            println("MINPACK Powell singular problem:")
        end
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        @test isapprox(x, [0.0, 0.0, 0.0, 0.0], atol=100*tol)

    end


    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
    @testset "... Test MINPACK Powell badly scaled problem" begin

        # Test function
        m = 2
        n = 2
        c1 = 1e4
        c2 = 1.0001
        function fx!(fvec, x)
            fvec[1] = c1*x[1]*x[2] - 1
            fvec[2] = exp(-x[1]) + exp(-x[2]) - c2
            return true
        end

        # Initial guess
        x = [0.0, 1.0]

        # Solver parameters
        tol = 1e-8
        maxiter = 50
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = false

        if loglevel >= info
            println("MINPACK Powell badly scaled problem:")
        end
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        @test isapprox(x, [1.0981593296997054e-5, 9.106146739867453], atol=100*tol)

    end


    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
    @testset "... Test MINPACK Wood problem" begin

        # Test function
        m = 4
        n = 4
        c3 = 2e2
        c4 = 2.02e1
        c5 = 1.98e1
        c6 = 1.8e2
        function fx!(F, x)
            temp1 = x[2] - x[1]^2
            temp2 = x[4] - x[3]^2
            F[1] = -c3*x[1]*temp1 - (1 - x[1])
            F[2] = c3*temp1 + c4*(x[2] - 1) + c5*(x[4] - 1)
            F[3] = -c6*x[3]*temp2 - (1 - x[3])
            F[4] = c6*temp2 + c4*(x[4] - 1) + c5*(x[2] - 1)
            return true
        end

        # Initial guess
        x = [-3.0, -1.0, -3.0, -1.0]

        # Solver parameters
        tol = 1e-8
        maxiter = 50
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = true

        if loglevel >= info
            println("MINPACK Wood problem:")
        end
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        @test isapprox(x, [-0.967974024937593, 0.9471391408178416, -0.9695163103315912, 0.9512476657923256], atol=100*tol)

    end


    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
    @testset "... Test MINPACK helical valley problem" begin

        # Test function
        m = 3
        n = 3
        tpi = 8*atan(1)
        c7 = 2.5e-1
        c8 = 5e-1
        function fx!(fvec, x)
            if x[1] > 0
                temp1 = atan(x[2]/x[1])/tpi
            elseif x[1] < 0
                temp1 = atan(x[2]/x[1])/tpi + c8
            else
                temp1 = c7*sign(x[2])
            end
            temp2 = sqrt(x[1]^2+x[2]^2)
            fvec[1] = 10(x[3] - 10*temp1)
            fvec[2] = 10(temp2 - 1)
            fvec[3] = x[3]
            return true
        end

        # Initial guess
        x = [-1.0, 0.0, 0.0]

        # Solver parameters
        tol = 1e-8
        maxiter = 50
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = false

        if loglevel >= info
            println("MINPACK helical valey problem:")
        end
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        @test isapprox(x, [1.0, 0.0, 0.0], atol=100*tol)

    end


    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
    @testset "... Test MINPACK Watson problem" begin

        # Solver parameters
        tol = 1e-7
        maxiter = 50
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = true

        for n in (6, 9)

            # Test function
            c9 = 2.9e1
            function fx!(fvec, x)
                fill!(fvec, 0)
                for i = 1:29
                    ti = i/c9
                    sum1 = 0.0
                    temp = 1.0
                    for j = 2:n
                        sum1 += (j-1)*temp*x[j]
                        temp *= ti
                    end
                    sum2 = 0.0
                    temp = 1.0
                    for j = 1:n
                        sum2 += temp*x[j]
                        temp *= ti
                    end
                    temp1 = sum1 - sum2^2 - 1
                    temp2 = 2*ti*sum2
                    temp = 1/ti
                    for k = 1:n
                        fvec[k] += temp*(k-1-temp2)*temp1
                        temp *= ti
                    end
                end
                temp = x[2] - x[1]^2 - 1
                fvec[1] += x[1]*(1-2temp)
                fvec[2] += temp
                return true
            end

            # Initial guess
            x = zeros(n)

            if loglevel >= info
                println("MINPACK Watson problem n = $n:")
            end
            y = zeros(Float64, n)
            fscale = ones(Float64, n)
            if loglevel == debug
                println("* start vector = $x")
                println("* scale vector = $fscale")
            end

            converged = solveNonlinearEquations!(fx!, n, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

            if loglevel == debug
                println("* solution = $x")
                fx!(y, x)
                println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
            end

            if n == 6
                @test isapprox(x, [-0.01572508640145857, 1.0124348693691099, -0.2329916259567373, 1.260430087799606, -1.5137289227222785, 0.9929964324311331], atol=100*tol)
            elseif n == 9
                @test isapprox(x, [-1.530703652140686e-5, 0.9997897039319482, 0.014763963693568192, 0.14634232829924446, 1.000821103005262, -2.617731140520362, 4.104403164480623, -3.143612278557626, 1.0526264080104843], atol=100*tol)
            end

        end

    end


#=
    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
    @testset "... Test MINPACK trigonometric problem" begin

        # Test function
        m = 10
        n = 10
        function fx!(fvec, x)
            for j = 1:n
                fvec[j] = cos(x[j])
            end
            sum1 = sum(fvec)
            for k = 1:n
                fvec[k] = n+k-sin(x[k]) - sum1 - k*fvec[k]
            end
            return true
        end

        # Start point
        x = ones(n)/n

        # Solver parameters
        tol = 1e-8
        maxiter = 500
        nonlinearity = extreme
        restricted = true
        quasi = false
        forwardjac = true

        if loglevel >= info
            println("MINPACK trigonometric problem:")
        end
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

    end
=#


    # Determined test from https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/test/minpack.jl
    @testset "... Test determined Rosenbrock problem" begin

        # Test function
        m = 2
        n = 2
        function fx!(fvec, x)
            fvec[1] = 1 - x[1]
            fvec[2] = 10*(x[2] - x[1]^2)
            return true
        end

        # Start point
        x = [0.0, -0.1]

        # Solver parameters
        tol = 1e-6
        maxiter = 500
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = false

        if loglevel >= info
            println("Determined Rosenbrock problem:")
        end
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        @test isapprox(x, [1.0, 1.0], atol=100*tol)

    end


    # Underdetermined test from https://de.wikipedia.org/wiki/Gau%C3%9F-Newton-Verfahren#Beispiel
    @testset "... Test underdetermined Rosenbrock problem" begin

        # Test function
        m = 1
        n = 2
        a = 1.0
        b = 100.0
        function fx!(fx::Vector, x::Vector)
            @assert(length(fx) == m)
            @assert(length(x) == n)
            fx[1] = (a - x[1])^2 + b*(x[2] - x[1]^2)^2
            return true
        end

        # Start point
        x = [0.0, -0.1]

        # Solver parameters
        tol = 1e-6
        maxiter = 500
        nonlinearity = high
        restricted = false
        quasi = false
        forwardjac = true

        if loglevel >= info
            println("Underdeterminded Rosenbrock problem:")
        end
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale")
        end

        converged = solveNonlinearEquations!(fx!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            fx!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        @test isapprox(x, [1.0, 1.0], atol=100*tol)

    end


    # Underdetermined slider-crank initial state problem 1
    @testset "... Test slider-crank initial state problem 1" begin

        # kinematic parameters
        wheelRadius = 0.1
        barLength = 0.5

        # function for computation of bar end y-component
        function res_z!(res::Vector, ang::Vector)
            phiWheel = ang[1]
            phiBar = ang[2]
            absPosJoint = wheelRadius*[-sin(phiWheel), cos(phiWheel)]
            relPosJointEnd = barLength*[cos(phiBar), sin(phiBar)]
            relTransJoint = [cos(phiWheel) -sin(phiWheel); sin(phiWheel) cos(phiWheel)]
            absPosEnd = absPosJoint + relTransJoint * relPosJointEnd
            res[1] = absPosEnd[2]
            return true
        end

        # initial guess
        phiWheel = 10.0
        phiBar = 5.0

        # solver parameters
        tol = 1e-6
        maxiter = 50
        log = true
        nonlinearity = mild
        restricted = false
        quasi = true
        forwardjac = false

        if loglevel >= info
            println("Slider-crank initial state problem 1:")
        end
        n = 2
        m = 1
        x = zeros(Float64, n)
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        x[1] = phiWheel
        x[2] = phiBar
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale\n")
        end

        converged = solveNonlinearEquations!(res_z!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            res_z!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        # consistent states
        phiWheel = x[1]
        phiBar   = x[2]

        @test isapprox([phiWheel, phiBar], [10.471751810654157, 5.1360050971720925], atol=100*tol)

    end


    # Underdetermined slider-crank initial state problem 2
    @testset "... Test slider-crank initial state problem 2" begin

        # kinematic parameters
        wheelRadius = 0.1
        barLength = 0.5

        # function for computation of joint position residual
        function res!(res::Vector, x::Vector)
            phiWheel = x[1]
            xBar = x[2]
            phiBar = x[3]
            absPosWheelJoint = wheelRadius*[-sin(phiWheel), cos(phiWheel)]
            absPosBarJoint = [xBar, 0.0] + barLength*[cos(phiBar), sin(phiBar)]
            res .= absPosWheelJoint .- absPosBarJoint
            return true
        end

        # initial guess
        phiWheel = 0.3
        xBar = -0.4
        phiBar = -0.1

        # solver parameters
        tol = 1e-6
        maxiter = 50
        log = true
        nonlinearity = high
        restricted = false
        quasi = true
        forwardjac = false

        if loglevel >= info
            println("Slider-crank initial state problem 2:")
        end
        n = 3
        m = 2
        x = zeros(Float64, n)
        y = zeros(Float64, m)
        fscale = ones(Float64, n)
        x[1] = phiWheel
        x[2] = xBar
        x[3] = phiBar
        if loglevel == debug
            println("* start vector = $x")
            println("* scale vector = $fscale\n")
        end

        converged = solveNonlinearEquations!(res!, m, x, fscale; nonlinearity=nonlinearity, restricted=restricted, xtol=tol, maxiter=maxiter, quasi=quasi, forwardjac=forwardjac, loglevel=loglevel)

        if loglevel == debug
            println("* solution = $x")
            res!(y, x)
            println("* rms(f(x*)) = $(norm(y)/sqrt(n))")
        end

        # consistent states
        phiWheel = x[1]
        xBar     = x[2]
        phiBar   = x[3]

        @test isapprox([phiWheel, xBar, phiBar], [0.30614037233377206, -0.5209621708411616, 0.1918759765115025], atol=100*tol)

    end

end


end
