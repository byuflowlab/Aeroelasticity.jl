"""
    Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

Two-dimensional aerodynamic model based on Wagner's function with state variables 
``\\lambda_1, \\lambda_2``, inputs ``u, v, \\omega``, and parameters ``a, b, a_0, 
\\alpha_0, c_{d0}, c_{m0}``
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

# --- Submodel Creation --- #

function Submodel(model::Wagner)

    C1 = model.C1
    C2 = model.C2
    eps1 = model.eps1
    eps2 = model.eps2

    # residual function
    fresid = (dx, x, y, p, t) -> wagner_residual(dx, x, y, p, t; C1, C2, eps1, eps2)

    # number of state, input, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(2)
    ny = Val(3)
    np = Val(6)

    # jacobian definitions
    ratejac = Identity()
    statejac = Linear((dx, x, y, p, t) -> wagner_state_jacobian(dx, x, y, p, t; model.eps1, model.eps2))
    inputjac = Nonlinear((dx, x, y, p, t) -> wagner_input_jacobian(dx, x, y, p, t; model.C1, model.C2, model.eps1, model.eps2))
    paramjac = Nonlinear((dx, x, y, p, t) -> wagner_parameter_jacobian(dx, x, y, p, t; model.C1, model.C2, model.eps1, model.eps2))
    tgrad = Zeros()

    # convenience functions for setting states, inputs, and parameters
    setstate = wagner_set_states!
    setinput = wagner_set_inputs!
    setparam = wagner_set_parameters!

    # convenience functions for separating states, inputs, and parameters
    sepstate = wagner_separate_states
    sepinput = wagner_separate_inputs
    sepparam = wagner_separate_parameters

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

# residual function
function wagner_residual(dx, x, y, p, t; C1, C2, eps1, eps2)
    # extract state rates
    d??1, d??2 = dx
    # extract states
    ??1, ??2 = x
    # extract inputs
    u, v, ?? = y
    # extract parameters
    a, b, a0, ??0, cd0, cm0 = p
    # return residual
    return SVector(
        d??1 + eps1*u/b*??1 - C1*eps1*u/b*(v + (b/2-a*b)*?? - u*??0),
        d??2 + eps2*u/b*??2 - C2*eps2*u/b*(v + (b/2-a*b)*?? - u*??0)
    )
end

# state jacobian function
function wagner_state_jacobian(dx, x, y, p, t; eps1, eps2)
    # extract inputs
    u, v, ?? = y
    # extract parameters
    a, b, a0, ??0, cd0, cm0 = p
    # return jacobian of residual with respect to the states
    return @SMatrix [eps1*u/b 0; 0 eps2*u/b]
end

# input jacobian function
function wagner_input_jacobian(dx, x, y, p, t; C1, C2, eps1, eps2)
    
    ??1, ??2 = x
    u, v, ?? = y
    a, b, a0, ??0, cd0, cm0 = p

    tmp1 = v/b + (1/2-a)*?? - 2*??0*u/b
    ??1dot_u = eps1/b*??1 - C1*eps1*tmp1
    ??2dot_u = eps2/b*??2 - C2*eps2*tmp1

    tmp2 = u/b
    ??1dot_v = -C1*eps1*tmp2
    ??2dot_v = -C2*eps2*tmp2

    tmp3 = u*(1/2-a)
    ??1dot_?? = -C1*eps1*tmp3
    ??2dot_?? = -C2*eps2*tmp3

    return @SMatrix [??1dot_u ??1dot_v ??1dot_??; ??2dot_u ??2dot_v ??2dot_??]
end

# parameter jacobian function
function wagner_parameter_jacobian(dx, x, y, p, t; C1, C2, eps1, eps2)

    ??1, ??2 = x
    u, v, ?? = y
    a, b, a0, ??0, cd0, cm0 = p

    ??1dot_a = C1*eps1*u*??
    ??1dot_b = -eps1*u/b^2*??1 + C1*eps1/b^2*(u*v - u^2*??0)
    ??1dot_a0 = 0
    ??1dot_??0 = C1*eps1*u^2/b
    ??1dot_cd0 = 0
    ??1dot_cm0 = 0

    ??2dot_a = C2*eps2*u*??
    ??2dot_b = -eps2*u/b^2*??2 + C2*eps2/b^2*(u*v - u^2*??0)
    ??2dot_a0 = 0
    ??2dot_??0 = C2*eps2*u^2/b
    ??2dot_cd0 = 0
    ??2dot_cm0 = 0

    return @SMatrix [
        ??1dot_a ??1dot_b ??1dot_a0 ??1dot_??0 ??1dot_cd0 ??1dot_cm0;
        ??2dot_a ??2dot_b ??2dot_a0 ??2dot_??0 ??2dot_cd0 ??2dot_cm0;
        ]
end

# convenience function for defining this model's state vector
wagner_set_states!(x; lambda) = x .= lambda

# convenience function for defining this model's input vector
function wagner_set_inputs!(y; u, v, omega)
    y[1] = u
    y[2] = v
    y[3] = omega
    return y
end

# convenience function for defining this model's parameter vector
function wagner_set_parameters!(p; a, b, a0, alpha0, cd0, cm0)
    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0
    p[5] = cd0
    p[6] = cm0
    return p
end

# convenience function for separating this model's state vector
wagner_separate_states(x) = (lambda = x,)

# convenience function for separating this model's input vector
wagner_separate_inputs(y) = (u=y[1], v=y[2], omega=y[3])

# convenience function for separating this model's parameter vector
wagner_separate_parameters(p) = (a=p[1], b=p[2], a0=p[3], alpha0=p[4], cd0=p[5], cm0=p[6])

# aerodynamic loads per unit span
function wagner_loads(a, b, ??, c, a0, ??0, cd0, cm0, C1, C2, u, v, ??, vdot, ??dot, ??1, ??2)
    # circulatory load factor
    tmp1 = a0*??*u*b
    # non-circulatory load factor
    tmp2 = pi*??*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ??0 = 1 - C1 - C2
    # Velocity Magnitude (squared)
    V2 = u^2 + v^2
    # Mach Number (squared)
    M2 = V2/c^2
    # Prandtl-Glauert correction factor
    beta = sqrt(1 - min(0.99, M2))
    # normal force at reference point
    N = tmp1*((v + d*?? - u*??0)*??0 + ??1 + ??2) + tmp2*(vdot/b + u/b*?? - a*??dot)
    # axial force at reference point
    A = -a0*??*b*((v + d*?? - u*??0)*??0 + ??1 + ??2)^2
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