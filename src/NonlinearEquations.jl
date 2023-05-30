"""
    module NonlinearEquations - Solve determined or underdetermined nonlinear equation systems

Initial implementation: Gerhard Hippmann, DLR-SR
"""
module NonlinearEquations

export mild, high, extreme, quiet, info, verbose, debug, solveNonlinearEquations!

using Printf
using LinearAlgebra
import ForwardDiff

@enum Nonlinearity begin
    mild
    high
    extreme
end

@enum Loglevel begin
    quiet
    info
    verbose
    debug
end

# Calculate RMS value (mean square root norm) of a vector
function rms(x::Vector)
    n = length(x)
    sum = 0
    for i = 1:n
       sum = sum + x[i]*x[i]
    end
    return sqrt(sum/n)
end


"""
    solveDeterminedNonlinearEquations!(F!::Function, x::Vector, scale::Union{Vector, Nothing}; nonlinearity::Nonlinearity = high, restricted = false, xtol = 1e-6, maxiter = 50, loglevel::Loglevel = info) -> convergence::Bool

Solve determined nonlinear equation system ``F(x) = 0`` by global Newton method with error oriented convergence criterion and adaptive trust region strategy.

`F!` is a C1 continuous function with identical length ``n`` of input and output. It has to be defined by
    F!(F::Vector, x::Vector) -> success::Bool
returning true in case of successful evaluation. Jacobian is computed by ForwardDiff package.

`x` is the vector of unknowns. On input it must contain an initial guess of the problem's solution, which is used as the start vector of the iteration. On output it contains the iteration result. If `solveNonlinearEquations!` returns true, the iteration converged and `x` contains an approximation of the solution which satisfies the error tolerance `xtol`.

`scale` is a vector of length ``n`` which contains positive scaling values used in computations of scaled norms and Jacobian scaling due to ``x_i / {scale}_i``. In case of `nothing` automatic scaling depends on `nonlinearity`.

`nonlinearity` defines the grade of non-linearity of the problem. Possible @enum values are `mild`, `high` and `extreme`. In case of `extreme` `restricted` is automatically set true.

`restricted` can be used to restrict the monotonicity test (i.e. damping factor) such that large iteration steps are avoided.

`xtol` is the error tolerance which the final approximate solution `x` must satisfy.

`maxiter` is the maximum number of allowed iterations.

'quasi' enables switch to local quasi Newton iteration (which avoids costly Jacobian evaluations) in case of advanced convergence.

`loglevel` is the log message level. Possible @enum values are `quiet`, `info`, `verbose` and `debug`.

This is a port of nleq_err.c available at http://elib.zib.de/pub/elib/codelib/NewtonLib/. Step reject is not implemented.

Reference: P. Deuflhard: Newton Methods for Nonlinear Problems. - Affine Invariance and Adaptive Algorithms. Series Computational Mathematics 35, Springer (2004)
"""
function solveDeterminedNonlinearEquations!( F!::Function,
                                             x::Vector,
                                             scale::Union{Vector, Nothing};
                                             nonlinearity::Nonlinearity = high,
                                             restricted = false,
                                             xtol = 1e-6,
                                             maxiter = 50,
                                             quasi = true,
                                             loglevel::Loglevel = info )

    n = length(x)
    @assert(n > 0)
    @assert(xtol > 0.0)
    @assert(maxiter > 0)

    # Constants
    SMALL = 1e-150
    LAMBDA_START_DEFAULT = 1e-2
    LAMBDA_MIN_DEFAULT = 1e-4
    LAMBDA_START_EXTREMELY_DEFAULT = 1e-4
    LAMBDA_MIN_EXTREMELY_DEFAULT = 1e-8
    THETA_MAX = 0.5

    # Work arrays
    xscale = similar(x, n)
    xthresh = similar(x, n)
    w = similar(x, n)
    dx = similar(x, n)
    dxbar = similar(x, n)
    fxk = similar(x, n)
    jac = similar(x, n, n)
    jcfg = ForwardDiff.JacobianConfig(F!, fxk, x)

    # Initialization
    k = 0
    cnvg = false
    nfcn = 0
    njac = 0
    normdx = 0.0
    normdxbar = 0.0
    precs = 0.0
    qnerr_iter = false

    # Set initial damping factor
    if nonlinearity == high
        lambda = LAMBDA_START_DEFAULT
        lambda_min = LAMBDA_MIN_DEFAULT
    elseif nonlinearity == extreme
        lambda = LAMBDA_START_EXTREMELY_DEFAULT
        lambda_min = LAMBDA_MIN_EXTREMELY_DEFAULT
        restricted = true
    else # mild
        lambda = 1.0
        lambda_min = LAMBDA_MIN_DEFAULT
    end

    # Set scaling vector
    if scale === nothing
        if nonlinearity == mild
            xscale .= 1.0
        else
            xscale .= xtol
        end
    else
        xscale .= scale
    end

    if loglevel >= verbose
        println("Solve determined nonlinear equations:")
        println("Problem dimension = $n")
        println("Prescribed relative precision = $xtol")
        println("The problem is specified as being $(nonlinearity)ly nonlinear.")
        println("The $(restricted ? "restricted" : "standard") monotonicity test will be applied.")
        println("The maximum permitted number of iteration steps is: $maxiter")
        if loglevel >= debug
            println("Start damping factor = $lambda")
            println("Start vector = $x")
            println("Scale vector = $xscale")
        end
    end

    # Evaluate F(x^0)
    okay = F!(fxk, x)
    nfcn = nfcn + 1
    if !okay
        error("Error in function evaluation $nfcn.")
    else
        normfk = rms(fxk)
        w .= x
        xthresh .= xscale

        if (loglevel >= debug)
            println("Start function vector = $fxk\n")
            println(" iter     norm_scl(dx)      norm(fk)    lambda")
        end
    end

    # For iteration index k = 0, 1,... do:
    while okay

        # Update scale vector
        #   xscale = max.(xthresh, max.(0.5*abs.(x) + abs.(w), SMALL))
        for i = 1:n
            xscale[i] = max(xthresh[i], max(0.5 * (abs(x[i]) + abs(w[i])), SMALL))
        end

        if k > 0
            # Recompute norms after rescaling
            normdxkm1 = rms(dx ./ xscale)
            normdxbar = rms(dxbar ./ xscale)
        end

        # 1. Step k: Evaluate Jacobian matrix F'(x^k)
        # In nleq_err.c fxk is not updated -> inconsistent after return from quasi Newton iteration
        ForwardDiff.jacobian!(jac, F!, fxk, x, jcfg)
        njac = njac + 1

        # Scale Jacobian and change sign of it
        for ic = 1:n
            jac[:,ic] = jac[:,ic] * -xscale[ic]
        end

        # Compute LU-factorization of Jacobian
        jaclu = lu(jac; check = false)
        if !issuccess(jaclu)
            println("LU factorization of Jacobian failed.")
            break
        end

        # Solve linear system -F'(x^k) * dx^k = F(x^k)
        ldiv!(dx, jaclu, fxk)

        # Descale predictor increment vector
        dx = dx .* xscale
        normdx = rms(dx ./ xscale)
        precs = normdx

        if loglevel >= debug
            @printf("%5i     %12e  %12e  %7f\n", k, normdx, normfk, lambda )
        end

        # Convergence test: If norm(dx^k) <= epsilon: Solution found: x^* := x^k + dx^k
        cnvg = normdx <= xtol
        if cnvg
            x .= x .+ dx
            k = k + 1
            break
        end

        # For k > 0: Compute a prediction value for the damping factor  (3.49)
        if k > 0
            w = dxbar .- dx
            s = rms(w ./ xscale)
            mue = (normdxkm1*normdxbar)/(s*normdx) * lambda
            lambda = min(1.0, mue)
        end
        reduced = false

        @label checkregularity
        # Regularity test: If lambda_k < lambda_min: Convergence failure
        if lambda < lambda_min
            if loglevel >= verbose
                println("Convergence failure, damping factor became too small.")
            end
            break
        end

        # 2. Compute the trial iterate x^{k+1} := x^k + lambda_k * dx^k  (3.43)
        xkp1 = x .+ lambda*dx

        # Evaluate F(x^{k+1})
        okay = F!(fxk, xkp1)
        nfcn = nfcn + 1
        if !okay
            println("Error in function evaluation $nfcn. Step reject not implemented.")
            break
        end
        normfkp1 = rms(fxk)

        # Solve linear system (old Jacobian, new right hand side) F'(x^k) * dxbar^{k+1} = -F(x^{k+1})  (page 138 4th formula)
        ldiv!(dxbar, jaclu, fxk)

        # Descale corrector increment vector
        dxbar = dxbar .* xscale
        normdxbar = rms(dxbar ./ xscale)
        precs = normdxbar

        if loglevel >= debug
            @printf("%5i  *  %12e  %12e  %7f\n", k, normdxbar, normfkp1, lambda )
        end

        # 3. Compute the monitoring quantities Theta_k := norm(dxbar^{k+1}) / norm(dx^k) and mu'_k := (0.5 * norm(dx^k) * lambda_k^2) / norm(dxbar_{k+1} - (1 - lambda_k) * dx^k)  (page 147 1st formula)
        theta = normdxbar / normdx  # Contraction factor  (page 148 1st formula)
        s = 1.0 - lambda
        w = dxbar .- s*dx
        mue = (0.5*normdx*lambda*lambda) / rms(w ./ xscale)

        # If Theta_k >= 1 (or, if restricted Theta_k > 1 - lambda_k/4)
        # then replace lambda_k by lambda'_k := min(mu'_k, lambda_k/2). Goto regularity test.
        # else let lambda'_k := min(1, lambda'_k)
        if (!restricted && theta >= 1.0) || (restricted && theta > 1.0-lambda/4.0)
            # Natural/restricted monotonicity test failed
            lambda_new = min(mue, 0.5*lambda)  # Corrected damping factor  (3.48)
            if lambda <= lambda_min
                lambda = lambda_new  # Does not make sense, bug in nleq_err.c?
            else
                lambda = max(lambda_new, lambda_min)
            end
            reduced = true
            @goto checkregularity
        else
            # Monotonicity test succeeded
            lambda_new = min(1.0, mue)  # (3.45)
        end

        # If lambda'_k == lambda_k == 1
        # then if norm(dxbar^{k+1}) <= epsilon: Solution found: x^* := x^{k+1} + dxbar^{k+1}
        #      else if Theta_k < 0.5 switch to quasi Newton method
        # else if lambda'_k >= 4*lambda_k replace lambda_k by lambda'_k and goto 2.
        # else accept x^{k+1} as new iterate: goto 1 with k := k + 1.
        if (lambda == 1.0) && (lambda_new == 1.0)
            # Iterates approach the solution
            cnvg = normdxbar <= xtol
            # Convergence test
            if cnvg
                x .= xkp1 .+ dxbar
                break
            end
            if quasi
                qnerr_iter = theta < THETA_MAX  # page 148 1st formula
            end
        elseif (lambda_new >= 4.0*lambda) && !reduced  # not documented?
            lambda = lambda_new
            @goto checkregularity
        end

        # Save previous iterate for scaling purposes and accept new iterate
        w .= x
        x .= xkp1
        normfk = normfkp1

        # Next step
        k = k + 1
        if k >= maxiter
            break
        end

        if qnerr_iter

            # Local error-oriented quasi Newton method using Broyden's 'good' rank-1 updates
            # Input  arguments: x, nfcn, normdx, jaclu, F!, n, scale, dx, maxiter
            # Output arguments: x, nfcn, normdx, normfk, cnvg, precs, okay, qnerr_iter
            # Overwrites: xscale

            # Work arrays
            dx_qn = similar(x, n, 1)
            sigma = similar(x, maxiter+2)
            fxk_qn = similar(x, n)
            v = similar(x, n)

            # Set scaling vector
            if scale === nothing
                xscale .= 1.0
            else
                xscale .= scale
            end

            # Initialization
            dx_qn[:,1] .= dx
            s = normdx
            sigma[1] = s * s

            # For k = 0, ..., k_max:
            k_qn = 0
            while true

                if k_qn > 0
                    # New iterate x_k+1
                    x .= x .+ dx_qn[:,k_qn+1]
                    precs = normdx
                    cnvg = sigma[k_qn+1] <= xtol * xtol
                    if cnvg
                        k_qn = k_qn + 1
                        break
                    end
                end

                # Allocate new column of dx_qn
                dx_qn = hcat(dx_qn, similar(x, n))

                # Evaluate F(x_k+1)
                okay = F!(fxk_qn, x)
                nfcn = nfcn + 1
                if !okay
                    println("Error in function evaluation $nfcn. Step reject not implemented.")
                    break
                end
                normfk = rms(fxk_qn)

                # Solve linear system J*v = -F_k+1
                ldiv!(v, jaclu, fxk_qn)

                # Descale v
                v = v .* xscale

                for i = 1:k_qn
                    alpha_bar = dot(v ./ xscale, dx_qn[:,i] ./ xscale)/n / sigma[i]
                    v = v .+ (alpha_bar * dx_qn[:,i+1])
                end

                alpha_kp1 = dot(v ./ xscale, dx_qn[:,k_qn+1] ./ xscale)/n / sigma[k_qn+1]
                thetak = sqrt(dot(v ./ xscale, v ./xscale)/n / sigma[k_qn+1])

                # Check contraction factor
                if thetak > THETA_MAX
                    if loglevel >= debug
                        @printf("%5i     %12e  %12e     QNERR  THETA!\n", k+k_qn, normdx, normfk )
                    end
                    k_qn = k_qn + 1
                    break
                end

                s = 1.0 - alpha_kp1
                dx_qn[:,k_qn+2] .= v / s
                sigma[k_qn+2] = dot(dx_qn[:,k_qn+2] ./ xscale, dx_qn[:,k_qn+2] ./ xscale)/n
                normdx = sqrt(sigma[k_qn+2])
                if loglevel >= debug
                    @printf("%5i     %12e  %12e     QNERR\n", k+k_qn, normdx, normfk )
                end
                k_qn = k_qn + 1

                if k + k_qn >= maxiter
                    break
                end

            end

            k = k + k_qn - 1  # Total number of iterations
            if cnvg
                break
            else
                normfk = normfkp1
                qnerr_iter = false
            end
            if k >= maxiter
                break
            end

        end

    end

    if k >= maxiter
        println("Max. allowed number of iterations = $maxiter exceeded.")
    end
    if loglevel >= info
        println("Solver $(cnvg ? "converged" : "did not converge").")
        if loglevel >= verbose
            println("Number of iterations = $k")
            println("Number of function evaluations = $nfcn")
            println("Number of Jacobian evaluations = $njac")
            if loglevel >= debug
               println("Solution vector = $x")
            end
            println("Precision = $precs")
        end
    end

    return cnvg

end


"""
    solveUnderdeterminedNonlinearEquations!(F!::Function, x::Vector, scale::Union{Vector, Nothing}; nonlinearity::Nonlinearity = high, restricted = false, xtol = 1e-6, maxiter = 50, quasi = true, loglevel::Loglevel = info) -> convergence::Bool

Solve nonlinear equation system ``F(x) = 0`` by global Gauss-Newton method with error oriented convergence criterion and adaptive trust region strategy.

`F!` is a C1 continuous function with length ``n`` of input and ``m`` of output vector where ``n >= m`` (determined or underdetermined system of equations). It has to be defined by
    F!(F::Vector, x::Vector) -> success::Bool
returning true in case of successful evaluation. Jacobian is computed by ForwardDiff package.

`m` is the length of the output vector of `F!`.

`x` is the vector of unknowns. On input it must contain an initial guess of the problem's solution, which is used as the start vector of the iteration. On output it contains the iteration result. If `solveUnderdeterminedNonlinearEquations!` returns true, the iteration converged and `x` contains an approximation of the solution which satisfies the error tolerance `xtol`.

`scale` is a vector of length ``n`` which contains positive scaling values used in computations of scaled norms and Jacobian scaling due to ``x_i / {scale}_i``. In case of `nothing` automatic scaling depends on `nonlinearity`.

`nonlinearity` defines the grade of non-linearity of the problem. Possible @enum values are `mild`, `high` and `extreme`. In case of `extreme` `restricted` is automatically set true.

`restricted` can be used to restrict the monotonicity test (i.e. damping factor) such that large iteration steps are avoided.

`xtol` is the error tolerance which the final approximate solution `x` must satisfy.

`maxiter` is the maximum number of allowed iterations.

'quasi' enables switch to local quasi Newton iteration (which avoids costly Jacobian evaluations) in case of advanced convergence.

`loglevel` is the log message level. Possible @enum values are `quiet`, `info`, `verbose` and `debug`.

Reference: P. Deuflhard: Newton Methods for Nonlinear Problems. - Affine Invariance and Adaptive Algorithms. Series Computational Mathematics 35, Springer (2004)
"""
function solveUnderdeterminedNonlinearEquations!( F!::Function,
                                                  m::Integer,
                                                  x::Vector,
                                                  scale::Union{Vector, Nothing};
                                                  nonlinearity::Nonlinearity = high,
                                                  restricted = false,
                                                  xtol = 1e-6,
                                                  maxiter = 50,
                                                  quasi = true,
                                                  loglevel::Loglevel = info )

    n = length(x)
    @assert(m > 0)
    @assert(n >= m)
    @assert(xtol > 0.0)
    @assert(maxiter > 0)

    # Constants
    SMALL = 1e-150
    LAMBDA_START_DEFAULT = 1e-2
    LAMBDA_MIN_DEFAULT = 1e-4
    LAMBDA_START_EXTREMELY_DEFAULT = 1e-4
    LAMBDA_MIN_EXTREMELY_DEFAULT = 1e-8
    THETA_MAX = 0.5
    DELKCOMP = false  # use [deltabar] = delta_k(x^{k-1}) instead of [deltabar] = 0

    # Work arrays
    xscale = similar(x, n)
    xthresh = similar(x, n)
    w = similar(x, n)
    dx = similar(x, n)
    dxbar = similar(x, n)
    fxk = similar(x, m)
    jac = similar(x, m, n)
    jcfg = ForwardDiff.JacobianConfig(F!, fxk, x)
    dxl = similar(x, n)
    Pdxbar = similar(x, n)
    fxl = similar(x, m)
    jacdxbar = similar(x, m)
    if DELKCOMP
        fxkm1 = similar(x, m)
        jacdxkm1 = similar(x, m)
    end

    # Initialization
    k = 0
    cnvg = false
    nfcn = 0
    njac = 0
    normdx = 0.0
    normdxbar = 0.0
    precs = 0.0
    qnerr_iter = false

    # Set initial damping factor
    if nonlinearity == high
        lambda = LAMBDA_START_DEFAULT
        lambda_min = LAMBDA_MIN_DEFAULT
    elseif nonlinearity == extreme
        lambda = LAMBDA_START_EXTREMELY_DEFAULT
        lambda_min = LAMBDA_MIN_EXTREMELY_DEFAULT
        restricted = true
    else # mild
        lambda = 1.0
        lambda_min = LAMBDA_MIN_DEFAULT
    end

    # Set scaling vector
    if scale === nothing
        if nonlinearity == mild
            xscale .= 1.0
        else
            xscale .= xtol
        end
    else
        xscale .= scale
    end

    if loglevel >= verbose
        println("Solve underdetermined nonlinear equations:")
        println("Problem dimension = $n")
        println("Prescribed relative precision = $xtol")
        println("The problem is specified as being $(nonlinearity)ly nonlinear.")
        println("The $(restricted ? "restricted" : "standard") monotonicity test will be applied.")
        println("The maximum permitted number of iteration steps is: $maxiter")
        if loglevel >= debug
            println("Start damping factor = $lambda")
            println("Start vector = $x")
            println("Scale vector = $xscale")
        end
    end

    # Evaluate F(x^0)
    okay = F!(fxk, x)
    nfcn = nfcn + 1
    if !okay
        error("Error in function evaluation $nfcn.")
    else
        normfk = rms(fxk)
        w .= x
        xthresh .= xscale

        if (loglevel >= debug)
            println("Start function vector = $fxk\n")
            println(" iter     norm_scl(dx)      norm(fk)    lambda")
        end
    end

    # For iteration index k = 0, 1,... do:
    while okay

        # Update scale vector
        #   xscale = max.(xthresh, max.(0.5*abs.(x) + abs.(w), SMALL))
        for i = 1:n
            xscale[i] = max(xthresh[i], max(0.5 * (abs(x[i]) + abs(w[i])), SMALL))
        end

        if k > 0
            # Recompute norms after rescaling
            normdxkm1 = rms(dx ./ xscale)
            normdxbar = rms(dxbar ./ xscale)
        end

        # 1. Step k: Evaluate Jacobian matrix F'(x^k)
        # In nleq_err.c fxk is not updated -> inconsistent after return from quasi Gauss-Newton iteration
        ForwardDiff.jacobian!(jac, F!, fxk, x, jcfg)
        njac = njac + 1

        # Scale Jacobian and change sign of it
        for ic = 1:n
            jac[:,ic] = jac[:,ic] * -xscale[ic]
        end

        # Compute QR-factorization of Jacobian
        jacqr = qr(jac, Val(true))

        # Solve linear system -F'(x^k) * dx^k = F(x^k)
        ldiv!(dx, jacqr, fxk)

        # Descale predictor increment vector
        dx = dx .* xscale
        normdx = rms(dx ./ xscale)
        precs = normdx

        if loglevel >= debug
            @printf("%5i     %12e  %12e  %7f\n", k, normdx, normfk, lambda )
        end

        # Convergence test: If norm(dx^k) <= epsilon: Solution found: x^* := x^k + dx^k
        cnvg = normdx <= xtol
        if cnvg
            x .= x .+ dx
            k = k + 1
            break
        end

        # For k > 0: Compute a prediction value for the damping factor
        if k > 0

            lambdakm1 = lambda  # Damping factor of last iterate

            # estimate for Lipschitz constant [omegabar_k]  (page 228 1st formula)
            # Solve -F'(x^k) * Pdxbar = -F'(x^k) * Delta(x^{k-1},x^k)  based on (4.79)
            mul!(jacdxbar, jac, dxbar ./ xscale)
            ldiv!(Pdxbar, jacqr, jacdxbar)  # Pdxbar = P(x^k) * Delta(x^{k-1},x^k)
            dxh = rms((dxbar .- dx) ./ xscale)^2 - rms((dxbar ./ xscale) .- Pdxbar)^2  # (Delta(x^{k-1},x^k) - Delta(x^k,x^k))^2 - ((I_n - P(x^k)) * Delta(x^{k-1},x^k))^2
            omegabar = sqrt(dxh) / (lambdakm1*normdxkm1*normdxbar)

            # Damping factor  (page 228 2nd formula)
            lambda = min(1.0, 1.0/(omegabar*normdx))

        end
        reduced = false

        @label checkregularity
        # Regularity test: If lambda_k < lambda_min: Convergence failure
        if lambda < lambda_min
            if loglevel >= verbose
                println("Convergence failure, damping factor became too small.")
            end
            break
        end

        # 2. Compute the trial iterate x^{k+1} := x^k + lambda_k * dx^k  (3.43)
        xkp1 = x .+ lambda*dx

        # Evaluate F(x^{k+1})
        okay = F!(fxk, xkp1)
        nfcn = nfcn + 1
        if !okay
            println("Error in function evaluation $nfcn. Step reject not implemented.")
            break
        end
        normfkp1 = rms(fxk)

        # Simplified Gauss-Newton correction
        # Solve linear system (old Jacobian, new right hand side) F'(x^k) * dxbar^{k+1} = F(x^{k+1})  page 199 3rd formula
        ldiv!(dxbar, jacqr, fxk)

        # Descale corrector increment vector
        dxbar = dxbar .* xscale
        normdxbar = rms(dxbar ./ xscale)
        precs = normdxbar

        if loglevel >= debug
            @printf("%5i  *  %12e  %12e  %7f\n", k, normdxbar, normfkp1, lambda )
        end

        # 3. Compute the monitoring quantity Theta_k := norm(dxbar^{k+1}) / norm(dx^k) and correction damping factor mue
        theta = normdxbar / normdx  # Contraction factor  (page 148 1st formula resp. page 213 1st formula resp. page 227 3rd formula)
        if k > 0
            if (lambdakm1 == 1.0) && (lambda == 1.0)
                # Iterates approach the solution
                del = theta  # page 215 1st paragraph
            elseif DELKCOMP
                # solve F'(x^k) * w = F(x^{k-1}) + F'(x^{k-1})*dx^{k-1}
                ldiv!(w, jacqr, fxkm1 + jacdxkm1)  # w = F'(x^k)^+ * r(x^{k-1})
                del1 = rms(w)
                # solve F'(x^k) * w = F(x^{k-1})
                ldiv!(w, jacqr, fxkm1)  # w = F'(x^k)^+ * F(x^{k-1})
                del2 = rms(w)
                # delta_k(x^{k-1})  (page 214 8th formula / page 217 1st formula)
                del = del1 / del2
            else
                del = 0.0
            end
        else
            del = 0.0
        end
        # Solve F'(x^k) *  dxl = F(x^k + lambda*Delta(x^k,x^k))
        w = x .+ lambda*dx
        okay = F!(fxl, w)
        nfcn = nfcn + 1
        if !okay
            println("Error in function evaluation $nfcn. Step reject not implemented.")
            break
        end
        ldiv!(dxl, jacqr, fxl)  # dxl = Delta(x^k,x^k+lambda*Delta(x^k,x^k))
        hk = 2.0*rms(dxl .- (1.0 - lambda)*dx ./ xscale) / (lambda*lambda*normdx)  # [h_k]  (4.68)
        mue = (1.0 - del) / hk  # Correction damping factor  (page 213 4th formula)

        # If Theta_k >= 1 (or, if restricted Theta_k > 1 - lambda_k/4)
        # then replace lambda_k by lambda'_k := min(mu'_k, lambda_k/2). Goto regularity test.
        # else let lambda'_k := min(1, lambda'_k)
        if (!restricted && theta >= 1.0) || (restricted && theta > 1.0-lambda/4.0)
            # Natural/restricted monotonicity test failed
            lambda_new = min(mue, 0.5*lambda)  # Corrected damping factor  (3.48)
            if lambda <= lambda_min
                lambda = lambda_new  # Does not make sense, bug in nleq_err.c?
            else
                lambda = max(lambda_new, lambda_min)
            end
            reduced = true
            @goto checkregularity
        else
            # Monotonicity test succeeded
            lambda_new = min(1.0, mue)  # (3.45)
        end

        # If lambda'_k == lambda_k == 1
        # then if norm(dxbar^{k+1}) <= epsilon: Solution found: x^* := x^{k+1} + dxbar^{k+1}
        #      else if Theta_k < 0.5 switch to quasi Gauss-Newton method
        # else if lambda'_k >= 4*lambda_k replace lambda_k by lambda'_k and goto 2.
        # else accept x^{k+1} as new iterate: goto 1 with k := k + 1.
        if (lambda == 1.0) && (lambda_new == 1.0)
            # Iterates approach the solution
            cnvg = normdxbar <= xtol
            # Convergence test
            if cnvg
                x .= xkp1 .+ dxbar
                break
            end
            if quasi
                qnerr_iter = theta < THETA_MAX  # page 148 1st formula
            end
        elseif (lambda_new >= 4.0*lambda) && !reduced  # not documented?
            lambda = lambda_new
            @goto checkregularity
        end

        # Save previous iterate for scaling purposes and accept new iterate
        w .= x
        x .= xkp1
        normfk = normfkp1
        if DELKCOMP
            fxkm1 .= fxk
            mul!(jacdxkm1, jac, dx ./ xscale)
        end

        # Next step
        k = k + 1
        if k >= maxiter
            break
        end

        if qnerr_iter

            # Local error-oriented quasi Gauss-Newton method using Broyden's 'good' rank-1 updates
            # Input  arguments: x, nfcn, normdx, jaclu, F!, n, scale, dx, maxiter
            # Output arguments: x, nfcn, normdx, normfk, cnvg, precs, okay, qnerr_iter
            # Overwrites: xscale

            # Work arrays
            dx_qn = similar(x, n, 1)
            sigma = similar(x, maxiter+2)
            fxk_qn = similar(x, m)
            v = similar(x, n)

            # Set scaling vector
            if scale === nothing
                xscale .= 1.0
            else
                xscale .= scale
            end

            # Initialization
            dx_qn[:,1] .= dx
            s = normdx
            sigma[1] = s * s

            # For k = 0, ..., k_max:
            k_qn = 0
            while true

                if k_qn > 0
                    # New iterate x_k+1
                    x .= x .+ dx_qn[:,k_qn+1]
                    precs = normdx
                    cnvg = sigma[k_qn+1] <= xtol * xtol
                    if cnvg
                        k_qn = k_qn + 1
                        break
                    end
                end

                # Allocate new column of dx_qn
                dx_qn = hcat(dx_qn, similar(x, n))

                # Evaluate F(x_k+1)
                okay = F!(fxk_qn, x)
                nfcn = nfcn + 1
                if !okay
                    println("Error in function evaluation $nfcn. Step reject not implemented.")
                    break
                end
                normfk = rms(fxk_qn)

                # Solve linear system J*v = -F_k+1
                ldiv!(v, jacqr, fxk_qn)

                # Descale v
                v = v .* xscale

                for i = 1:k_qn
                    alpha_bar = dot(v ./ xscale, dx_qn[:,i] ./ xscale)/n / sigma[i]
                    v = v .+ (alpha_bar * dx_qn[:,i+1])
                end

                alpha_kp1 = dot(v ./ xscale, dx_qn[:,k_qn+1] ./ xscale)/n / sigma[k_qn+1]
                s = 1.0 - alpha_kp1
                dx_qn[:,k_qn+2] .= v / s
                thetak = rms(dx_qn[:,k_qn+2]) / rms(dx_qn[:,k_qn+1])  # page 225 3rd formula

                # Check contraction factor
                if thetak > THETA_MAX
                    if loglevel >= debug
                        @printf("%5i     %12e  %12e     QNERR  THETA!\n", k+k_qn, normdx, normfk )
                    end
                    k_qn = k_qn + 1
                    break
                end

                sigma[k_qn+2] = dot(dx_qn[:,k_qn+2] ./ xscale, dx_qn[:,k_qn+2] ./ xscale)/n
                normdx = sqrt(sigma[k_qn+2])
                if loglevel >= debug
                    @printf("%5i     %12e  %12e     QNERR\n", k+k_qn, normdx, normfk )
                end
                k_qn = k_qn + 1

                if k + k_qn >= maxiter
                    break
                end

            end

            k = k + k_qn - 1  # Total number of iterations
            if cnvg
                break
            else
                normfk = normfkp1
                qnerr_iter = false
            end
            if k >= maxiter
                break
            end

        end

    end

    if k >= maxiter
        println("Max. allowed number of iterations = $maxiter exceeded.")
    end
    if loglevel >= info
        println("Solver $(cnvg ? "converged" : "did not converge").")
        if loglevel >= verbose
            println("Number of iterations = $k")
            println("Number of function evaluations = $nfcn")
            println("Number of Jacobian evaluations = $njac")
            if loglevel >= debug
               println("Solution vector = $x")
            end
            println("Precision = $precs")
        end
    end

    return cnvg

end


"""
    solveNonlinearEquations!(F!::Function,
                             m::Integer,
                             x::Vector,
                             scale::Union{Vector, Nothing};
                             nonlinearity::Nonlinearity = high,
                             restricted = false,
                             xtol = 1e-6,
                             maxiter = 50,
                             quasi = true,
                             loglevel::Loglevel = info) -> convergence::Bool

Solve nonlinear equation system ``F(x) = 0`` by global Newton or Gauss-Newton method with error oriented convergence criterion and adaptive trust region strategy. Optionally Broyden's 'good' Jacobian rank-1 updates are used.

`F!` is a C1 continuous function with length ``n`` of input and ``m`` of output vector where ``n >= m`` (determined or underdetermined system of equations). It has to be defined by
    F!(F::Vector, x::Vector) -> success::Bool
returning true in case of successful evaluation. Jacobian is computed by [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package.

`m` is the length of the output vector of `F!`.

`x` is the vector of unknowns. On input it must contain an initial guess of the problem's solution, which is used as the start vector of the iteration. On output it contains the iteration result. If `solveNonlinearEquations!` returns true, the iteration converged and `x` contains an approximation of the solution which satisfies the error tolerance `xtol`.

`scale` is a vector of length ``n`` which contains positive scaling values used in computations of scaled norms and Jacobian scaling due to ``x_i / {scale}_i``. In case of `nothing` automatic scaling depends on `nonlinearity`.

`nonlinearity` defines the grade of non-linearity of the problem. Possible @enum values are `mild`, `high` and `extreme`. In case of `extreme` `restricted` is automatically set true.

`restricted` can be used to restrict the monotonicity test (i.e. damping factor) such that large iteration steps are avoided.

`xtol` is the error tolerance which the final approximate solution `x` must satisfy.

`maxiter` is the maximum number of allowed iterations.

`quasi` enables switch to local quasi Newton iteration (which avoids costly Jacobian evaluations) in case of advanced convergence.

`loglevel` is the log message level. Possible @enum values are `quiet`, `info`, `verbose` and `debug`.

Reference: P. Deuflhard: Newton Methods for Nonlinear Problems. - Affine Invariance and Adaptive Algorithms. Series Computational Mathematics 35, Springer (2004)

In case of ``n = m``, i.e. determined system of equations, the Newton method described in textbook section 3.3.3 (port of [NLEQ_ERR](https://www.zib.de/codelib/NewtonLib/)) is used. In case of ``n > m``, i.e. underdetermined system of equations, the Gauss-Newton method described in textbook section 4.4 is used to compute the 'shortest' least square solution of ``F(x) = 0``.
"""
function solveNonlinearEquations!(F!::Function,
                                  m::Integer,
                                  x::Vector,
                                  scale::Union{Vector, Nothing};
                                  nonlinearity::Nonlinearity = high,
                                  restricted = false,
                                  xtol = 1e-6,
                                  maxiter = 50,
                                  quasi = true,
                                  loglevel::Loglevel = info)

    n = length(x)
    @assert(m > 0)
    @assert(n >= m)
    @assert(xtol > 0.0)
    @assert(maxiter > 0)

    if n == m
        if (loglevel >= debug)
            println("Call solver for determined nonlinear equation system with n = m = $n")
        end
        solveDeterminedNonlinearEquations!(F!, x, scale; nonlinearity=nonlinearity, restricted=restricted, xtol=xtol, maxiter=maxiter, quasi=quasi, loglevel=loglevel)
    else
        if (loglevel >= debug)
            println("Call solver for underdetermined nonlinear equation system with n = $n and m = $m")
        end
        solveUnderdeterminedNonlinearEquations!(F!, m, x, scale; nonlinearity=nonlinearity, restricted=restricted, xtol=xtol, maxiter=maxiter, quasi=quasi, loglevel=loglevel)
    end

end

end
