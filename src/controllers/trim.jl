"""
    Trim{N} <: AbstractModel

Trimmed aircraft controller model with state variables corresponding to control
inputs which are used to trim the aircraft and inputs corresponding to the total
forces and moments on the aircraft ``F_x, F_y, F_z, M_x, M_y, M_z``.
"""
struct Trim{N} <: AbstractModel
    state_indices::NTuple{N,Int}
    force_indices::NTuple{N,Int}
end

"""
    Trim(state_indices::NTuple{N}, force_indices::NTuple{N})

Initialize an object of type [`Trim`](@ref). Use the state variables
corresponding to `state_indices` to trim the forces/moments corresponding to
`force_indices`.
"""
Trim(state_indices, force_indices)

Trim() = Trim((1,2,3,4,5,6), (1,2,3,4,5,6))

# --- Traits --- #

number_of_states(::Type{Trim{N}}) where N = N
number_of_inputs(::Type{Trim{N}}) where N = N
number_of_parameters(::Type{<:Trim}) = 0

inplaceness(::Type{<:Trim}) = OutOfPlace()

rate_jacobian_type(::Type{<:Trim}) = Zeros()
state_jacobian_type(::Type{<:Trim}) = Zeros()
input_jacobian_type(::Type{<:Trim}) = Identity()
parameter_jacobian_type(::Type{<:Trim}) = Empty()
time_gradient_type(::Type{<:Trim}) = Zeros()

# --- Methods --- #

get_residual(::Trim, dx, x, y, p, t) = y

# --- Performance Overloads --- #

# NOTE: No performance overloads are needed since the jacobians are fully
# defined by their provided types

# --- Convenience Methods --- #

set_states!(x, model::Trim; delta) = x .= delta

function set_inputs!(y, model::Trim; Fx, Fy, Fz, Mx, My, Mz)

    y[1] = Fx
    y[2] = Fy
    y[3] = Fz
    y[4] = Mx
    y[5] = My
    y[6] = Mz

    return y
end

set_parameters!(p, model::Trim) = p

separate_states(model::Trim, x) = (delta = x,)

function separate_inputs(model::Trim, y)

    return (Fx = y[1], Fy = y[2], Fz = y[3], Mx = y[4], My = y[5], Mz = y[6])
end

separate_parameters(model::Trim, p) = ()

# --- Internal Methods for Couplings --- #
