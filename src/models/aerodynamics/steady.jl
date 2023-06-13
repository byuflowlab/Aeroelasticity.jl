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
function steady_loads(a, b, a0, α0, cd0, cm0, rho, beta, u, v)
    # normal force at reference point
    N = a0*rho*b*u*(v - u*α0)
    # axial force at reference point
    A = -a0*rho*b*(v - u*α0)^2
    # moment at reference point
    M = 2*rho*b^2*u^2*cm0 + (b/2 + a*b)*N
    # apply compressibility correction (use constant beta to preserve linearity)
    N = N / beta
    A = A / beta
    M = M / beta
    # add skin friction drag (note that Prandtl-Glauert correction applies only to pressure)
    A += rho*b*u^2*cd0
    return SVector(N, A, M)
end