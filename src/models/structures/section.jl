"""
    Section()

Construct a typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
function Section()

    # residual function
    fresid = section_residual

    # number of state, input, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(4)
    ny = Val(2)
    np = Val(5)

    # jacobian definitions
    ratejac = Constant(section_rate_jacobian)
    statejac = Constant(section_state_jacobian)
    inputjac = Invariant(section_input_jacobian)
    paramjac = Linear(section_parameter_jacobian)
    tgrad = Zeros()

    # convenience functions for setting states, inputs, and parameters
    setstate = section_setstate!
    setinput = section_setinput!
    setparam = section_setparam!

    # convenience functions for separating states, inputs, and parameters
    sepstate = section_sepstate
    sepinput = section_sepinput
    sepparam = section_sepparam

    # model definition
    return Model{false}(fresid, nx, ny, np;
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

# residual function
function section_residual(dx, x, y, p, t)

    dh, dθ, dhdot, dθdot = dx
    h, θ, hdot, θdot = x
    L, M = y
    kh, kθ, m, Sθ, Iθ = p   

    r1 = dh - hdot
    r2 = dθ - θdot
    r3 = m*dhdot + Sθ*dθdot + kh*h + L
    r4 = Sθ*dhdot + Iθ*dθdot + kθ*θ - M

    return SVector(r1, r2, r3, r4)
end

# rate jacobian function
function section_rate_jacobian(p)
    kh, kθ, m, Sθ, Iθ = p
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m Sθ; 0 0 Sθ Iθ]
end

# state jacobian function
function section_state_jacobian(p)
    kh, kθ, m, Sθ, Iθ = p
    return @SMatrix [0 0 -1 0; 0 0 0 -1; kh 0 0 0; 0 kθ 0 0]
end

# input jacobian definition
const section_input_jacobian = @SMatrix [0 0; 0 0; 1 0; 0 -1]

# parameter jacobian function
function section_parameter_jacobian(dx, x, y, p, t)
    dh, dθ, dhdot, dθdot = dx
    h, θ, hdot, θdot = x
    return @SMatrix [0 0 0 0 0; 0 0 0 0 0; h 0 dhdot dθdot 0; 0 θ 0 dhdot dθdot]
end

# convenience function for defining this model's state vector
function section_setstate!(x; h, theta, hdot, thetadot)
    x[1] = h
    x[2] = theta
    x[3] = hdot
    x[4] = thetadot
    return x
end

# convenience function for defining this model's input vector
function section_setinput!(y; L, M)
    y[1] = L
    y[2] = M
    return y
end

# convenience function for defining this model's parameter vector
function section_setparam!(p; kh, ktheta, m, Stheta, Itheta)
    p[1] = kh
    p[2] = ktheta
    p[3] = m
    p[4] = Stheta
    p[5] = Itheta
    return p
end

# convenience function for separating this model's state vector
section_sepstate(x) = (h = x[1], theta = x[2], hdot = x[3], thetadot = x[4])

# convenience function for separating this model's input vector
section_sepinput(y) = (L = y[1], M = y[2])

# convenience function for separating this model's parameter vector
section_sepparam(p) = (kh = p[1], ktheta = p[2], m = p[3], Stheta = p[4], Itheta = p[5])

# --- Internal Methods for Couplings with this Model --- #

# airfoil local linear/angular velocities
function section_velocities(U, θ, hdot, θdot)

    u = U
    v = U*θ + hdot
    ω = θdot

    return SVector(u, v, ω)
end

# airfoil local linear/angular velocities derivatives
section_velocities_U(θ) = SVector(1, θ, 0)
section_velocities_θ(U) = SVector(0, U, 0)
section_velocities_hdot() = SVector(0, 1, 0)
section_velocities_θdot() = SVector(0, 0, 1)

# airfoil local linear/angular accelerations
function section_accelerations(dhdot, dθdot)
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    return SVector(udot, vdot, ωdot)
end

# airfoil local linear/angular acceleration derivatives
section_accelerations_dhdot() = SVector(0, 1, 0)
section_accelerations_dθdot() = SVector(0, 0, 1)
