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

number_of_states(::Type{<:LiftingLineSection}) = 6
number_of_inputs(::Type{<:LiftingLineSection}) = 6
number_of_parameters(::Type{<:LiftingLineSection}) = 0
inplaceness(::Type{LiftingLineSection}) = OutOfPlace()
