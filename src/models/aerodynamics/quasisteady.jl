"""
    QuasiSteady

Two-dimensional aerodynamic model based on quasi-steady thin airfoil theory with parameters 
``a, b, a_0, \\alpha_0, c_{d0}, c_{m0}``.
"""
struct QuasiSteady end

"""
    QuasiSteady()

Initialize a model of type [`QuasiSteady`](@ref)
"""
QuasiSteady()

# --- Submodel Creation --- #

function Submodel(::QuasiSteady)

    # number of parameters (use Val(N) to use inferrable dimensions)
    np = Val(6)

    # convenience function for setting parameters
    setparam = quasisteady_set_parameters!

    # convenience function for separating parameters
    sepparam = quasisteady_separate_parameters

    # model definition
    return Submodel(np; setparam, sepparam)
end

# --- Internal Methods --- #

# convenience function for defining this model's parameter vector
function quasisteady_set_parameters!(p; a, b, a0, alpha0, cd0, cm0)
    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0
    p[5] = cd0
    p[6] = cm0
    return p
end

# convenience function for separating this model's parameter vector
quasisteady_separate_parameters(p) = (a=p[1], b=p[2], a0=p[3], alpha0=p[4], cd0=p[5], cm0=p[6])

# aerodynamic loads per unit span
function quasisteady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Velocity Magnitude (squared)
    V2 = u^2 + v^2
    # Mach Number (squared)
    M2 = V2/c^2
    # Prandtl-Glauert correction factor
    beta = sqrt(1 - min(0.99, M2))
    # normal force at reference point
    N = tmp1*(v + d*ω - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # axial force at reference point
    A = -a0*ρ*b*(v + d*ω - u*α0)^2
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + (b/8 - a*b/2)*ωdot) + 2*ρ*b^2*u^2*cm0 + (b/2 + a*b)*N
    # apply compressibility correction
    N = N / beta
    A = A / beta
    M = M / beta
    # add skin friction drag
    A += ρ*b*u^2*cd0
    # return loads
    return SVector(N, A, M)
end