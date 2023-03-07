"""
    Steady

Two-dimensional aerodynamic model based on steady thin airfoil theory with parameters 
``a, b, a_0, \\alpha_0, c_{d0}, c_{m0}``.
"""
struct Steady end

"""
    Steady()

Initialize a model of type [`Steady`](@ref)
"""
Steady()

# residual function (no states so it is satisfied by default)
(::Steady)(resid, dx, x, y, p, t) = resid .= 0

number_of_states(::Steady) = 0

number_of_parameters(::Steady) = 6

# --- Internal Methods --- #

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