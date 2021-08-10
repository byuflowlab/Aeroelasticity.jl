"""
    LiftingLineSectionControl <: AbstractModel

Lifting line section control input model with state variable ``\\delta``, zero
inputs, and zero parameters.  Two-dimensional control surface models may be
extended to three dimensional models by coupling with this model.  Note that this
model has no rate equations of its own since its state variables are defined as
functions of the 3D system's control inputs.
"""
struct LiftingLineSectionControl <: AbstractModel end

number_of_states(::Type{<:LiftingLineSectionControl}) = 1
number_of_inputs(::Type{<:LiftingLineSectionControl}) = 0
number_of_parameters(::Type{<:LiftingLineSectionControl}) = 0
inplaceness(::Type{<:LiftingLineSectionControl}) = OutOfPlace()
