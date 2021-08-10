"""
    Trim{N,TI} <: AbstractModel

Trimmed aircraft controller model with up to six state variables and inputs.
State variables for this model correspond to the control inputs which are used
to trim the aircraft.  Inputs for this model correspond to the subset of the
total forces and moments ``F_x, F_y, F_z, M_x, M_y, M_z`` which are trimmed by
this model.
"""
struct Trim{N,TI} <: AbstractModel
    state_indices::NTuple{N,TI}
    force_indices::NTuple{N,TI}
end

"""
    Trim(state_indices::NTuple{N,TI}, force_indices::NTuple{N,TI})

Initialize an object of type [`Trim`](@ref). Use the state variables
corresponding to `state_indices` to trim the forces/moments corresponding to
`force_indices`.
"""
Trim(state_indices, force_indices)

Trim() = Trim((1,2,3,4,5,6), (1,2,3,4,5,6))

# --- Traits --- #

number_of_states(::Type{Trim{N,TI}}) where {N,TI} = N
number_of_inputs(::Type{Trim{N,TI}}) where {N,TI} = N
number_of_parameters(::Type{<:Trim}) = 0
inplaceness(::Type{<:Trim}) = OutOfPlace()
mass_matrix_type(::Type{<:Trim}) = Zeros()
state_jacobian_type(::Type{<:Trim}) = Zeros()
input_jacobian_type(::Type{<:Trim}) = Identity()

# --- Methods --- #

get_rates(::Trim, x, y, p, t) = y

# --- Unit Testing Methods --- #

get_lhs(::Trim, dx, x, y, p, t) = zero(y)
