"""
    Peters(N)

Construct an aerodynamic model based on Peters' finite state model with `N` state variables, 
inputs ``u, \\omega, \\dot{v}, \\dot{\\omega}`` and parameters ``a, b, a_0, \\alpha_0, 
c_{d,0}, c_{m,0}``
"""
function Peters(N)

    # model constants (which may be used when this model is coupled with other models)
    Abar, bbar, cbar = peters_constants(N)

    # residual function
    fresid = (dx, x, y, p, t) -> peters_residual(dx, x, y, p, t; Abar, cbar)

    # number of state, input, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(N)
    ny = Val(4)
    np = Val(6)

    # jacobian definitions
    ratejac = Invariant(Abar)
    statejac = Linear((dx, x, y, p, t) -> peters_state_jacobian(dx, x, y, p, t; cbar))
    inputjac = Nonlinear((dx, x, y, p, t) -> peters_input_jacobian(dx, x, y, p, t; Abar, cbar))
    paramjac = Nonlinear((dx, x, y, p, t) -> peters_parameter_jacobian(dx, x, y, p, t; Abar, cbar))
    tgrad = Zeros()

    # convenience functions for setting states, inputs, and parameters
    setstate = peters_setstate!
    setinput = peters_setinput!
    setparam = peters_setparam!

    # convenience functions for separating states, inputs, and parameters
    sepstate = peters_sepstate
    sepinput = peters_sepinput
    sepparam = peters_sepparam

    # model definition
    return Model{false}(fresid, nx, ny, np;
        constants = (Abar=Abar, bbar=bbar, cbar=cbar),
        ratejac = ratejac,
        statejac = statejac,
        inputjac = inputjac,
        paramjac = paramjac,
        tgrad = tgrad,
        setstate = setstate,
        setinput = setinput,
        setparam = setparam,
        sepstate = sepstate,
        sepinput = sepinput,
        sepparam = sepparam,
    )
end

# --- Internal Methods for this Model --- #

# function for defining model constants
function peters_constants(N)

    b = zeros(N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(big(N + n - 1))/factorial(big(N - n - 1))*
            1/factorial(big(n))^2
    end
    b[N] = (-1)^(N-1)

    c = zeros(N)
    for n = 1:N
        c[n] = 2/n
    end

    d = zeros(N)
    d[1] = 1/2

    D = zeros(N, N)
    for m in 1:N-1
        n = m + 1
        D[n, m] = 1/(2*n)
    end
    for m in 2:N
        n = m - 1
        D[n, m] = -1/(2*n)
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return SMatrix{N,N}(A), SVector{N}(b), SVector{N}(c)
end

# residual function
function peters_residual(dx, x, y, p, t; Abar, cbar)

    dλ = similar_type(cbar, eltype(dx))(dx)
    λ = similar_type(cbar, eltype(x))(x)
    u, ω, vdot, ωdot = y
    a, b, a0, α0, cd0, cm0 = p
    
    return Abar*dλ - cbar*(vdot + u*ω + (b/2-a*b)*ωdot) + u/b*λ
end

# rate jacobian function
peters_rate_jacobian(; Abar) = Abar

# state jacobian function
function peters_state_jacobian(dx, x, y, p, t; cbar)

    u, ω, vdot, ωdot = y
    a, b, a0, α0, cd0, cm0 = p

    return u/b*Diagonal(ones(similar_type(cbar)))
end

# input jacobian function
function peters_input_jacobian(dx, x, y, p, t; Abar, cbar)

    λ = similar_type(cbar, eltype(x))(x)
    u, ω, vdot, ωdot = y
    a, b, a0, α0, cd0, cm0 = p

    dλ_u = -cbar*ω + λ/b
    dλ_ω = -u*cbar
    dλ_vdot = -cbar
    dλ_ωdot = -(b/2-a*b)*cbar

    return hcat(dλ_u, dλ_ω, dλ_vdot, dλ_ωdot)
end

# parameter jacobian function
function peters_parameter_jacobian(dx, x, y, p, t; Abar, cbar)

    λ = similar_type(cbar, eltype(x))(x)
    u, ω, vdot, ωdot = y
    a, b, a0, α0, cd0, cm0 = p

    dλ_a = b*ωdot*cbar
    dλ_b = -cbar*(1/2-a)*ωdot - u/b^2*λ
    dλ_a0 = zero(cbar)
    dλ_α0 = zero(cbar)
    dλ_cd0 = zero(cbar)
    dλ_cm0 = zero(cbar)

    return hcat(dλ_a, dλ_b, dλ_a0, dλ_α0, dλ_cd0, dλ_cm0)
end

# convenience function for defining this model's state vector
peters_setstate!(x; lambda) = x .= lambda

# convenience function for defining this model's input vector
function peters_setinput!(y; u, omega, vdot, omegadot)
    y[1] = u
    y[2] = omega
    y[3] = vdot
    y[4] = omegadot
    return y
end

# convenience function for defining this model's parameter vector
function peters_setparam!(p; a, b, a0, alpha0, cd0, cm0)
    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0
    p[5] = cd0
    p[6] = cm0
    return p
end

# convenience function for separating this model's state vector
peters_sepstate = (x) -> (lambda=x,)

# convenience function for separating this model's input vector
peters_sepinput = (y) -> (u=y[1], omega=y[2], vdot=y[3], omegadot=y[4])

# convenience function for separating this model's parameter vector
peters_sepparam = (p) -> (a=p[1], b=p[2], a0=p[3], alpha0=p[4], cd0=p[5], cm0=p[6])

# --- Internal Methods for Couplings with this Model --- #

# aerodynamic loads per unit span
function peters_loads(a, b, ρ, c, a0, α0, cd0, cm0, bbar, u, v, ω, vdot, ωdot, λ)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # Velocity Magnitude (squared)
    V2 = u^2 + v^2
    # Mach Number (squared)
    M2 = V2/c^2
    # Prandtl-Glauert correction factor
    beta = sqrt(1 - min(0.99, M2))
    # normal force at the reference point
    N = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # axial force at the reference point
    A = -a0*ρ*b*(v + d*ω - λ0 - u*α0)^2
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