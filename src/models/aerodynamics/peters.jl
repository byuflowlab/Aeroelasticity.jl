"""
    Peters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}}

Two-dimensional aerodynamic model based on Peters' finite state model with `N` state 
variables, inputs ``u, \\omega, \\dot{v}, \\dot{\\omega}`` and parameters ``a, b, a_0, 
\\alpha_0, c_{d,0}, c_{m,0}``
"""
struct Peters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}}
    A::TA
    b::TV
    c::TV
end

"""
    Peters{N, TF=Float64}()

Initialize a model with type [`Peters`](@ref) with `N` state variables
"""
Peters()

Peters{N}() where N = Peters{N,Float64}()

function Peters{N,TF}() where {N,TF}

    A, b, c = peters_constants(N, TF)

    return Peters(A, b, c)
end

# --- Submodel Creation --- #

function Submodel(model::Peters{N}) where N

    Abar = model.A
    bbar = model.b
    cbar = model.c

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
    setstate = peters_set_states!
    setinput = peters_set_inputs!
    setparam = peters_set_parameters!

    # convenience functions for separating states, inputs, and parameters
    sepstate = peters_separate_states
    sepinput = peters_separate_inputs
    sepparam = peters_separate_parameters

    # model definition
    return Submodel{false}(fresid, nx, ny, np;
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

# --- Internal Methods --- #

# function for defining model constants
function peters_constants(N, TF)

    b = zeros(TF, N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(big(N + n - 1))/factorial(big(N - n - 1))*
            1/factorial(big(n))^2
    end
    b[N] = (-1)^(N-1)

    c = zeros(TF, N)
    for n = 1:N
        c[n] = 2/n
    end

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
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

    d?? = similar_type(cbar, eltype(dx))(dx)
    ?? = similar_type(cbar, eltype(x))(x)
    u, ??, vdot, ??dot = y
    a, b, a0, ??0, cd0, cm0 = p
    
    return Abar*d?? - cbar*(vdot + u*?? + (b/2-a*b)*??dot) + u/b*??
end

# rate jacobian function
peters_rate_jacobian(; Abar) = Abar

# state jacobian function
function peters_state_jacobian(dx, x, y, p, t; cbar)

    u, ??, vdot, ??dot = y
    a, b, a0, ??0, cd0, cm0 = p

    return u/b*Diagonal(ones(similar_type(cbar)))
end

# input jacobian function
function peters_input_jacobian(dx, x, y, p, t; Abar, cbar)

    ?? = similar_type(cbar, eltype(x))(x)
    u, ??, vdot, ??dot = y
    a, b, a0, ??0, cd0, cm0 = p

    d??_u = -cbar*?? + ??/b
    d??_?? = -u*cbar
    d??_vdot = -cbar
    d??_??dot = -(b/2-a*b)*cbar

    return hcat(d??_u, d??_??, d??_vdot, d??_??dot)
end

# parameter jacobian function
function peters_parameter_jacobian(dx, x, y, p, t; Abar, cbar)

    ?? = similar_type(cbar, eltype(x))(x)
    u, ??, vdot, ??dot = y
    a, b, a0, ??0, cd0, cm0 = p

    d??_a = b*??dot*cbar
    d??_b = -cbar*(1/2-a)*??dot - u/b^2*??
    d??_a0 = zero(cbar)
    d??_??0 = zero(cbar)
    d??_cd0 = zero(cbar)
    d??_cm0 = zero(cbar)

    return hcat(d??_a, d??_b, d??_a0, d??_??0, d??_cd0, d??_cm0)
end

# convenience function for defining this model's state vector
peters_set_states!(x; lambda) = x .= lambda

# convenience function for defining this model's input vector
function peters_set_inputs!(y; u, omega, vdot, omegadot)
    y[1] = u
    y[2] = omega
    y[3] = vdot
    y[4] = omegadot
    return y
end

# convenience function for defining this model's parameter vector
function peters_set_parameters!(p; a, b, a0, alpha0, cd0, cm0)
    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0
    p[5] = cd0
    p[6] = cm0
    return p
end

# convenience function for separating this model's state vector
peters_separate_states = (x) -> (lambda=x,)

# convenience function for separating this model's input vector
peters_separate_inputs = (y) -> (u=y[1], omega=y[2], vdot=y[3], omegadot=y[4])

# convenience function for separating this model's parameter vector
peters_separate_parameters = (p) -> (a=p[1], b=p[2], a0=p[3], alpha0=p[4], cd0=p[5], cm0=p[6])

# aerodynamic loads per unit span
function peters_loads(a, b, ??, c, a0, ??0, cd0, cm0, bbar, u, v, ??, vdot, ??dot, ??)
    # circulatory load factor
    tmp1 = a0*??*u*b
    # non-circulatory load factor
    tmp2 = pi*??*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    ??0 = 1/2 * bbar'*??
    # Velocity Magnitude (squared)
    V2 = u^2 + v^2
    # Mach Number (squared)
    M2 = V2/c^2
    # Prandtl-Glauert correction factor
    beta = sqrt(1 - min(0.99, M2))
    # normal force at the reference point
    N = tmp1*(v + d*?? - ??0 - u*??0) + tmp2*(vdot/b + u/b*?? - a*??dot)
    # axial force at the reference point
    A = -a0*??*b*(v + d*?? - ??0 - u*??0)^2
    # moment at reference point
    M = -tmp2*(vdot/2 + u*?? + b*(1/8 - a/2)*??dot) + 2*??*b^2*u^2*cm0 + (b/2 + a*b)*N
    # apply compressibility correction
    N = N / beta
    A = A / beta
    M = M / beta
    # add skin friction drag
    A += ??*b*u^2*cd0
    # return loads
    return SVector(N, A, M)
end