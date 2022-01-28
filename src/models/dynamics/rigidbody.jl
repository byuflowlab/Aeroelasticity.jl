"""
    rigidbody_model(state_indices=(), rate_indices=(), prescribed_values=())

Construct a six-degree of freedom rigid-body model with state variables ``x, y, z, \\phi,
\\theta, \\psi, u, v, w, p, q, r``, inputs ``m, I_{xx}, I_{yy}, I_{zz}, I_{xz},
I_{xy}, I_{yz}, F_x, F_y, F_z, M_x, M_y, M_z``, and zero parameters.  Use the state 
variables corresponding to `state_indices` to set the state rates corresponding to
`rate_indices` to the values specified in `prescribed_values`. Otherwise, allow the state 
variables and their respective rates to be defined by their rate equations.
"""
function rigidbody_model(state_indices=(), rate_indices=(), prescribed_values=())
    
    # residual function
    fresid = (dx, x, y, p, t) -> rigidbody_residual(dx, x, y, p, t; 
        state_indices, rate_indices, prescribed_values)

    # number of state, input, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(12)
    ny = Val(13)
    np = Val(0)

    # jacobian definitions
    ratejac = iszero(length(state_indices)) ? Identity() : Linear() # TODO: define ratejac function
    statejac = Nonlinear() # TODO: define statejac function
    inputjac = Nonlinear() # TODO: define inputjac function
    paramjac = Nonlinear() # TODO: define paramjac function
    tgrad = Zeros()

    # convenience functions for setting states, inputs, and parameters
    setstate = rigidbody_setstate!
    setinput = rigidbody_setinput!
    setparam = rigidbody_setparam!

    # convenience functions for separating states, inputs, and parameters
    sepstate = rigidbody_sepstate
    sepinput = rigidbody_sepinput
    sepparam = rigidbody_sepparam
    
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
function rigidbody_residual(dstates, states, inputs, parameters, time;
    state_indices, rate_indices, prescribed_values)
    
    # extract states
    x, y, z, ϕ, θ, ψ, u, v, w, p, q, r = states
    
    # extract inputs
    m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz = inputs
    
    # rigid body kinematics
    xdot, ydot, zdot, ϕdot, θdot, ψdot = rigidbody_kinematics(x, y, z, ϕ, θ, ψ, u, v, w, p, q, r)
    
    # rigid body dynamics
    udot, vdot, wdot, pdot, qdot, rdot = rigidbody_dynamics(u, v, w, p, q, r, m,
        Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz)
    
    # calculated state rates
    rates = SVector(xdot, ydot, zdot, ϕdot, θdot, ψdot, udot, vdot, wdot, pdot, qdot, rdot)
    
    # default residual
    residual = dstates - rates
    
    # replace elements of the residual vector with new constraints
    for (istate, irate, value) in zip(state_indices, rate_indices, prescribed_values)
        # replace state rate equation for `istate` with the prescribed constraint
        residual = setindex(residual, value - rates[irate], istate)
    end
    
    # return result
    return residual
end

# rigid body kinematics
function rigidbody_kinematics(x, y, z, ϕ, θ, ψ, u, v, w, p, q, r)

    Vb = SVector(u, v, w)
    Ωb = SVector(p, q, r)

    sϕ, cϕ = sincos(ϕ)
    sθ, cθ = sincos(θ)
    sψ, cψ = sincos(ψ)

    # linear kinematics
    Rib = @SMatrix [cθ*cψ    cθ*sψ           -sθ;
         sϕ*sθ*cψ - cϕ*sψ sϕ*sθ*sψ + cϕ*cψ sϕ*cθ;
         cϕ*sθ*cψ + sϕ*sψ cϕ*sθ*sψ - sϕ*cψ cϕ*cθ]
    rdot = Rib' * Vb

    # angular kinematics
    ϕdot = p + (q*sϕ + r*cϕ)*sθ/cθ
    θdot = q*cϕ - r*sϕ
    ψdot = (q*sϕ + r*cϕ)/cθ

    return SVector(rdot..., ϕdot, θdot, ψdot)
end

# rigid body dynamics
function rigidbody_dynamics(u, v, w, p, q, r, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz,
    Fx, Fy, Fz, Mx, My, Mz)

    F = SVector(Fx, Fy, Fz)
    M = SVector(Mx, My, Mz)

    Vb = SVector(u, v, w)
    Ωb = SVector(p, q, r)

    # linear dynamics
    vdot = F/m - cross(Ωb, Vb)

    # angular dynamics
    I = @SMatrix [Ixx -Ixy -Ixz; -Iyz Iyy -Iyz; -Ixz -Iyz Izz]
    ωdot = I \ (M - cross(Ωb, I*Ωb))

    return SVector(vdot..., ωdot...)
end

# convenience function for defining this model's state vector
function rigidbody_setstate!(states; x, y, z, ϕ, θ, ψ, u, v, w, p, q, r)

    states[1] = x
    states[2] = y
    states[3] = z
    states[4] = ϕ
    states[5] = θ
    states[6] = ψ
    states[7] = u
    states[8] = v
    states[9] = w
    states[10] = p
    states[11] = q
    states[12] = r

    return states
end

# convenience function for defining this model's input vector
function rigidbody_setinput!(inputs; m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz)

    inputs[1] = m
    inputs[2] = Ixx
    inputs[3] = Iyy
    inputs[4] = Izz
    inputs[5] = Ixz
    inputs[6] = Ixy
    inputs[7] = Iyz
    inputs[8] = Fx
    inputs[9] = Fy
    inputs[10] = Fz
    inputs[11] = Mx
    inputs[12] = My
    inputs[13] = Mz

    return inputs
end

# convenience function for defining this model's parameter vector
rigidbody_setparam!(p) = p

# convenience function for separating this model's state vector
rigidbody_sepstate(states) = (x = states[1], y = states[2], z = states[3],
    ϕ = states[4], θ = states[5], ψ = states[6],
    u = states[7], v = states[8], w = states[9],
    p = states[10], q = states[11], r = states[12])

# convenience function for separating this model's input vector
rigidbody_sepinput(inputs) = (m = inputs[1],
    Ixx = inputs[2], Iyy = inputs[3], Izz = inputs[4],
    Ixz = inputs[5], Ixy = inputs[6], Iyz = inputs[7],
    Fx = inputs[8], Fy = inputs[9], Fz = inputs[10],
    Mx = inputs[11], My = inputs[12], Mz = inputs[13])

# convenience function for separating this model's parameter vector
rigidbody_sepparam(parameters) = ()

# --- Internal Methods for Couplings with this Model --- #

function rigidbody_velocities(rigidbody_states)

    V = SVector(rigidbody_states.u, rigidbody_states.v, rigidbody_states.w)
    Omega = SVector(rigidbody_states.p, rigidbody_states.q, rigidbody_states.r)
    
    return V, Omega
end

function rigidbody_accelerations(rigidbody_rates)

    a = SVector(rigidbody_rates.u, rigidbody_rates.v, rigidbody_rates.w)
    alpha = SVector(rigidbody_rates.p, rigidbody_rates.q, rigidbody_rates.r)

    return a, alpha
end