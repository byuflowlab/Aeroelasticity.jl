"""
    Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

Two-dimensional aerodynamic model based on Wagner's function with state variables 
``\\lambda_1, \\lambda_2``, inputs ``u, v, \\omega``, and parameters ``a, b, a_0, 
\\alpha_0, c_{d0}, c_{m0}``.
"""
struct Wagner{TF}
    C1::TF
    C2::TF
    eps1::TF
    eps2::TF
end

"""
    Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

Initialize a model of type [`Wagner`](@ref)
"""
Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3) = Wagner(promote(C1, C2, eps1, eps2)...)

# residual function
function (wagner::Wagner)(resid, dx, x, y, p, t)
    # extract constants
    @unpack C1, C2, eps1, eps2 = wagner
    # extract state rates
    dλ1, dλ2 = dx
    # extract states
    λ1, λ2 = x
    # extract inputs
    u, v, ω = y
    # extract parameters
    a, b, a0, α0, cd0, cm0 = p
    # define residual
    resid[1] = dλ1 + eps1*u/b*λ1 - C1*eps1*u/b*(v + (b/2-a*b)*ω - u*α0)
    resid[2] = dλ2 + eps2*u/b*λ2 - C2*eps2*u/b*(v + (b/2-a*b)*ω - u*α0)
    # return residual
    return resid
end

number_of_states(::Wagner) = 2

number_of_parameters(::Wagner) = 6

# --- Internal Methods --- #

# aerodynamic loads per unit span
function wagner_loads(a, b, ρ, c, a0, α0, cd0, cm0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # Velocity Magnitude (squared)
    V2 = u^2 + v^2
    # Mach Number (squared)
    M2 = V2/c^2
    # Prandtl-Glauert correction factor
    beta = sqrt(1 - min(0.99, M2))
    # normal force at reference point
    N = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # axial force at reference point
    A = -a0*ρ*b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)^2
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + 2*ρ*b^2*u^2*cm0 + (b/2 + a*b)*N
    # apply compressibility correction
    N = N / beta
    A = A / beta
    M = M / beta
    # add skin friction drag
    A += ρ*b*u^2*cd0
    # return loads
    return SVector(N, A, M)
end