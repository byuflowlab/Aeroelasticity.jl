"""
    DynamicStall{TI,TV} <: AbstractModel

Dynamic stall model constructed by coupling an invisicid flow model and a
viscous flow model
"""
struct DynamicStall{TI,TV} <: AbstractModel
    inviscid::TI # inviscid (attached flow) model
    viscous::TV # viscous flow model
end

# Viscous Flow Models
struct Onera{TF}
    a0::TF = 0.3
    a2::TF = 0.2
    e2::TF = -2.86
    r0::TF = 0.2
    r2::TF = 0.2
end

# --- Constructors --- #

"""
    DynamicStall(inviscid, viscous)

Initialize an object of type [`DynamicStall`](@ref)
"""
DynamicStall()

# --- Traits --- #

function number_of_states(::Type{DynamicStall{TI,TV}}) where {TI,TV}
    return number_of_states(TI) + number_of_states(TV)
end
function number_of_inputs(::Type{DynamicStall{TI,TV}}) where {TI,TV}
    return number_of_inputs(TI) + number_of_inputs(TV)
end
function number_of_parameters(::Type{DynamicStall{TI,TV}}) where {TI,TV}
    return number_of_parameters(TI) + number_of_parameters(TV)
end
inplaceness(::Type{<:DynamicStall}) = OutOfPlace()
mass_matrix_type(::Type{<:DynamicStall}) = Identity()
state_jacobian_type(::Type{<:DynamicStall}) = Linear()
input_jacobian_type(::Type{<:DynamicStall}) = Nonlinear()

# --- Methods --- #

function get_rates(model::Wagner, λ, d, p, t) where {N,TF,SV,SA}
    # extract states
    λ1, λ2 = λ
    # extract inputs
    u, v, ω = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # calculate rates
    return wagner_rates(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)
end

# --- Performance Overloads --- #

function get_state_jacobian(model::Wagner, λ, d, p, t)
    # extract inputs
    u, v, ω = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    ε1 = model.eps1
    ε2 = model.eps2
    # jacobian with respect to aerodynamic states
    return wagner_state_jacobian(a, b, ε1, ε2, u)
end

function get_input_jacobian(model::Wagner, λ, d, p, t)
    # extract states
    λ1, λ2 = λ
    # extract inputs
    u, v, ω = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # return jacobian
    return wagner_input_jacobian(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)
end

# --- Unit Testing Methods --- #

get_lhs(::Wagner, dλ, λ, d, p, t) = dλ

# --- Convenience Methods --- #

set_states!(y, model::Wagner; lambda) = y .= lambda

function set_inputs!(y, model::Wagner; u, v, omega)

    y[1] = u
    y[2] = v
    y[3] = omega

    return y
end

function set_parameters!(p, model::Wagner; a, b, a0, alpha0)

    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0

    return p
end

function separate_states(model::Wagner, x)

    return (lambda = x,)
end

function separate_inputs(model::Wagner, y)

    return (u=y[1], v=y[2], omega=y[3])
end

function separate_parameters(model::Wagner, p)

    return (a=p[1], b=p[2], a0=p[3], alpha0=p[4])
end

# --- Internal Methods --- #
