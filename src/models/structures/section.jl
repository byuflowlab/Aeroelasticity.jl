"""
    Section()

Typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
struct Section end

# --- Submodel Creation --- #

function Submodel(::Section)

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
    setstate = section_set_states!
    setinput = section_set_inputs!
    setparam = section_set_parameters!

    # convenience functions for separating states, inputs, and parameters
    sepstate = section_separate_states
    sepinput = section_separate_inputs
    sepparam = section_separate_parameters

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

# --- Additional Functions --- #

function section_coordinates(h, θ; 
    a=0, 
    b=0.5,
    xcoord = [1.0, 0.993844, 0.975528, 0.945503, 0.904508, 0.853553, 0.793893, 0.726995, 
        0.654508, 0.578217, 0.5, 0.421783, 0.345492, 0.273005, 0.206107, 0.146447, 
        0.095492, 0.054497, 0.024472, 0.006156, 0.0, 0.006156, 0.024472, 0.054497, 
        0.095492, 0.146447, 0.206107, 0.273005, 0.345492, 0.421783, 0.5, 0.578217, 
        0.654508, 0.726995, 0.793893, 0.853553, 0.904508, 0.945503, 0.975528, 0.993844, 1.0],
    ycoord = [0.00126, 0.00212, 0.004642, 0.008658, 0.013914, 0.020107, 0.026905, 0.033962, 
        0.040917, 0.047383, 0.05294, 0.057148, 0.059575, 0.059848, 0.057714, 0.053083, 
        0.046049, 0.036867, 0.025893, 0.013503, 0.0, -0.013503, -0.025893, -0.036867, 
        -0.046049, -0.053083, -0.057714, -0.059848, -0.059575, -0.057148, -0.05294, 
        -0.047383, -0.040917, -0.033962, -0.026905, -0.020107, -0.013914, -0.008658, 
        -0.004642, -0.00212, -0.00126],
    )

    xplot = similar(xcoord)
    yplot = similar(ycoord)
    for i = 1:length(xcoord)
        xplot[i] = (xcoord[i] - 0.5 - a)*2*b*cos(θ) - ycoord[i]*2*b*sin(θ)
        yplot[i] = (xcoord[i] - 0.5 - a)*2*b*sin(θ) + ycoord[i]*2*b*cos(θ) + h
    end

    return xplot, yplot
end

# --- Internal Methods --- #

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
function section_set_states!(x; h, theta, hdot, thetadot)
    x[1] = h
    x[2] = theta
    x[3] = hdot
    x[4] = thetadot
    return x
end

# convenience function for defining this model's input vector
function section_set_inputs!(y; L, M)
    y[1] = L
    y[2] = M
    return y
end

# convenience function for defining this model's parameter vector
function section_set_parameters!(p; kh, ktheta, m, Stheta, Itheta)
    p[1] = kh
    p[2] = ktheta
    p[3] = m
    p[4] = Stheta
    p[5] = Itheta
    return p
end

# convenience function for separating this model's state vector
section_separate_states(x) = (h = x[1], theta = x[2], hdot = x[3], thetadot = x[4])

# convenience function for separating this model's input vector
section_separate_inputs(y) = (L = y[1], M = y[2])

# convenience function for separating this model's parameter vector
section_separate_parameters(p) = (kh = p[1], ktheta = p[2], m = p[3], Stheta = p[4], Itheta = p[5])

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
