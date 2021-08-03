"""
    Trim{Fx,Fy,Fz,Mx,My,Mz} <: AbstractModel

Controller which modifies control surface deflections to maintain trimmed
operating conditions with state variables ``\delta_1, \delta_2, \delta_3,
\delta_4, \delta_5, \delta_6`` and inputs ``F_x, F_y, F_z, M_x, M_y, M_z``.  A
subset of these state variables and inputs may also be used.
"""
struct Trim{Fx,Fy,Fz,Mx,My,Mz} <: AbstractModel end

"""
    Trim()

Initialize an object of type [`Trim`](@ref)
"""
Trim() = Trim{true, true, true, true, true, true}

"""
    Trim{Fx,Fy,Fz,Mx,My,Mz}()

Initialize an object of type [`Trim`](@ref)
"""
Trim{Fx,Fy,Fz,Mx,My,Mz}()

# --- Traits --- #

function number_of_states(::Type{Trim{Fx,Fy,Fz,Mx,My,Mz}}) where {Fx,Fy,Fz,Mx,My,Mz}

    return count((Fx, Fy, Fz, Mx, My, Mz))
end

function number_of_inputs(::Type{Trim{Fx,Fy,Fz,Mx,My,Mz}}) where {Fx,Fy,Fz,Mx,My,Mz}

    return count((Fx, Fy, Fz, Mx, My, Mz))
end

number_of_parameters(::Type{<:Trim}) = 0

inplaceness(::Type{<:Trim}) = OutOfPlace()

mass_matrix_type(::Type{<:Trim}) = Zeros()

state_jacobian_type(::Type{<:Trim}) = Zeros()

input_jacobian_type(::Type{<:Trim}) = Identity()

# --- Methods --- #

get_rates(::Trim, x, y, p, t) = y

# --- Unit Testing Methods --- #

get_lhs(::Trim, dx, x, y, p, t) = zero(y)
