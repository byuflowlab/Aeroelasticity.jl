"""
    LiftingLineSection

Lifting line section model with state variables ``v_x, v_y, v_z, \\omega_x, \\omega_y, 
\\omega_z``, inputs ``f_x, f_y, f_z, m_x, m_y, m_z``, and parameters ``\\rho, c``.  

Note that this model doesn't define any residuals and is defined purely as a convenience
to allow aerodynamic models to be more easily used with [`LiftingLine`](@ref).  Also note
that in order to preserve compatibility with the default coupling functions for the
[`LiftingLine`](@ref) model, no additional parameters may be introduced when creating a
coupling function involving this model.
"""
struct LiftingLineSection end

"""
    LiftingLineSection()

Initialize a model of type [`LiftingLineSection`](@ref)
"""
LiftingLineSection()

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
