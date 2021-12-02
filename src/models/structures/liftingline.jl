"""
    LiftingLineSection()

Construct a lifting line section model with state variables ``v_x, v_y, v_z, \\omega_x,
\\omega_y, \\omega_z``, inputs ``F_x', F_y', F_z', M_x, M_y, M_z``, and parameters 
``\\rho, c``.  Two-dimensional aerodynamic models may be extended to three-dimensional 
models by coupling with this model.  Note that this model has no rate equations of its own 
since its state variables are defined as functions of the 3D structural model's state 
variables.  Also note that additional parameters may not be used when coupling with this 
model, this is a limitation of the [`LiftingLine`](@ref) model.
"""
function LiftingLineSection()
   
    # residual function (no residual function)
    f = nothing 

    # number of state, input, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(6)
    ny = Val(6)
    np = Val(2)

    # convenience functions for setting states, inputs, and parameters
    setstate = liftinglinesection_setstate!
    setinput = liftinglinesection_setinput!
    setparam = liftinglinesection_setparam!

    # convenience functions for separating states, inputs, and parameters
    sepstate = liftinglinesection_sepstate
    sepinput = liftinglinesection_sepinput
    sepparam = liftinglinesection_sepparam

    return Model{false}(f, nx, ny, np)
end

# --- Internal Methods for this Model --- #

# convenience function for defining this model's state vector
function liftinglinesection_setstate!(x; v, omega)
    x[1] = v[1]
    x[2] = v[2]
    x[3] = v[3]
    x[4] = omega[1]
    x[5] = omega[2]
    x[6] = omega[3]
    return x
end

# convenience function for defining this model's input vector
function liftinglinesection_setinput!(y; F, M)
    y[1] = F[1]
    y[2] = F[2]
    y[3] = F[3]
    y[4] = M[1]
    y[5] = M[2]
    y[6] = M[3]
    return y
end

# convenience function for defining this model's parameter vector
function liftinglinesection_setparam!(p; rho, c)
    p[1] = rho
    p[2] = c
    return p
end

# convenience function for separating this model's state vector
liftinglinesection_sepstate(x) = (v=view(x, 1:3), omega=view(x, 4:6))

# convenience function for separating this model's input vector
liftinglinesection_sepinput(y) = (F=view(y, 1:3), M=view(y, 4:6))

# convenience function for separating this model's parameter vector
liftinglinesection_sepparam(p) = (rho=p[1], c=p[2])

# --- Internal Methods for Couplings with this Model --- #

# local freestream linear/angular velocities
function liftinglinesection_velocities(vx, vz, ωy)
    u = vx
    v = vz
    ω = ωy
    return u, v, ω
end

# local freestream linear/angular accelerations
function liftinglinesection_accelerations(dvx, dvz, dωy)
    udot = dvx
    vdot = dvz
    ωdot = dωy
    return udot, vdot, ωdot
end
