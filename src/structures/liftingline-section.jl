"""
    LiftingLineSection <: AbstractModel

Lifting line section model with state variables ``v_x, v_y, v_z, \\omega_x,
\\omega_y, \\omega_z``, inputs ``F_x', F_y', F_z', M_x, M_y, M_z``, and zero
parameters.  Two-dimensional aerodynamic models may be extended to
three-dimensional models by coupling with this model.  Note that
this model has no rate equations of its own since its state variables are
defined as functions of the 3D structural model's state variables.
"""
struct LiftingLineSection <: AbstractModel end

"""
    LiftingLineSection()

Initialize an object of type [`LiftingLineSection`](@ref)
"""
LiftingLineSection()

# --- Traits --- #

number_of_states(::Type{<:LiftingLineSection}) = 6
number_of_inputs(::Type{<:LiftingLineSection}) = 6
number_of_parameters(::Type{<:LiftingLineSection}) = 0

inplaceness(::Type{LiftingLineSection}) = OutOfPlace()

# --- Internal Methods for Couplings --- #

function liftingline_velocities(vx, vz, ωy)
    u = vx
    v = vz
    ω = ωy
    return u, v, ω
end

function liftingline_accelerations(dvx, dvz, dωy)
    udot = dvx
    vdot = dvz
    ωdot = dωy
    return udot, vdot, ωdot
end
