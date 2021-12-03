"""
    steady_model()

Construct a 2D aerodynamic model based on steady thin airfoil theory with parameters 
``a, b, a_0, \\alpha_0, c_{d0}, c_{m0}``.
"""
function steady_model()

    # number of parameters (use Val(N) to use inferrable dimensions)
    np = Val(6)

    # convenience function for setting parameters
    setparam = steady_setparam!

    # convenience function for separating parameters
    sepparam = steady_sepparam

    # model definition
    return NoStateModel(np; setparam, sepparam)
end

# --- Internal Methods for this Model --- #

# convenience function for defining this model's parameter vector
function steady_setparam!(p; a, b, a0, alpha0, cd0, cm0)
    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0
    p[5] = cd0
    p[6] = cm0
    return p
end

# convenience function for separating this model's parameter vector
steady_sepparam(p) = (a=p[1], b=p[2], a0=p[3], alpha0=p[4], cd0=p[5], cm0=p[6])

# --- Internal Methods for Couplings with this Model --- #

# aerodynamic loads per unit span
function steady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v)
    # Velocity Magnitude (squared)
    V2 = u^2 + v^2
    # Mach Number (squared)
    M2 = V2/c^2
    # Prandtl-Glauert correction factor
    beta = sqrt(1 - min(0.99, M2))
    # normal force at reference point
    N = a0*ρ*b*u*(v - u*α0)
    # axial force at reference point
    A = -a0*ρ*b*(v - u*α0)^2
    # moment at reference point
    M = 2*ρ*b^2*u^2*cm0 + (b/2 + a*b)*N
    # apply compressibility correction
    N = N / beta
    A = A / beta
    M = M / beta
    # add skin friction drag
    A += ρ*b*u^2*cd0
    return SVector(N, A, M)
end