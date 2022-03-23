"""
    LiftingLineSection

Lifting line section model with state variables ``v_x, v_y, v_z, \\omega_x,
\\omega_y, \\omega_z``, inputs ``f_x, f_y, f_z, m_x, m_y, m_z``, and parameters 
``\\rho, c``.  Two-dimensional aerodynamic models may be extended to three-dimensional 
models by coupling with this model.  Note that this model has no rate equations of its own 
since its state variables are defined as functions of the 3D structural model's state 
variables.  Also note that additional parameters may not be used when coupling with this 
model, this limitation is enforced in order to ensure compatability with the 
[`LiftingLine`](@ref) model.
"""
struct LiftingLineSection end

"""
    LiftingLineSection()

Initialize a model of type [`LiftingLineSection`](@ref)
"""
LiftingLineSection()

# --- Submodel Creation --- #

function Submodel(::LiftingLineSection)
   
    # residual function (no residual function)
    fresid = nothing

    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(6)
    ny = Val(6)
    np = Val(2)

    # convenience functions for setting states, inputs, and parameters
    setstate = liftingline_section_setstate!
    setinput = liftingline_section_setinput!
    setparam = liftingline_section_setparam!

    # convenience functions for separating states, inputs, and parameters
    sepstate = liftingline_section_sepstate
    sepinput = liftingline_section_sepinput
    sepparam = liftingline_section_sepparam

    return Submodel{false}(fresid, nx, ny, np;
        setstate = setstate,
        setinput = setinput,
        setparam = setparam,
        sepstate = sepstate,
        sepinput = sepinput,
        sepparam = sepparam)
end

# --- Internal Methods for this Model --- #

# convenience function for defining this model's state vector
function liftingline_section_setstate!(x; v, omega)
    x[1] = v[1]
    x[2] = v[2]
    x[3] = v[3]
    x[4] = omega[1]
    x[5] = omega[2]
    x[6] = omega[3]
    return x
end

# convenience function for defining this model's input vector
function liftingline_section_setinput!(y; f, m)
    y[1] = f[1]
    y[2] = f[2]
    y[3] = f[3]
    y[4] = m[1]
    y[5] = m[2]
    y[6] = m[3]
    return y
end

# convenience function for defining this model's parameter vector
function liftingline_section_setparam!(p; rho, c)
    p[1] = rho
    p[2] = c
    return p
end

# convenience function for separating this model's state vector
liftingline_section_sepstate(x) = (v=view(x, 1:3), omega=view(x, 4:6))

# convenience function for separating this model's input vector
liftingline_section_sepinput(y) = (f=view(y, 1:3), m=view(y, 4:6))

# convenience function for separating this model's parameter vector
liftingline_section_sepparam(p) = (rho=p[1], c=p[2])

# --- Internal Methods for Couplings with this Model --- #

# local freestream linear/angular velocities
function liftingline_section_velocities(vx, vz, ωy)
    u = vx
    v = vz
    ω = ωy
    return u, v, ω
end

# local freestream linear/angular accelerations
function liftingline_section_accelerations(dvx, dvz, dωy)
    udot = dvx
    vdot = dvz
    ωdot = dωy
    return udot, vdot, ωdot
end
