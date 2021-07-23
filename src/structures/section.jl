"""
    TypicalSection <: AbstractModel

Typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
struct TypicalSection <: AbstractModel end

"""
    TypicalSection()

Initialize an object of type [`TypicalSection`](@ref)
"""
TypicalSection()

# --- Traits --- #

number_of_states(::Type{TypicalSection}) = 4
number_of_inputs(::Type{TypicalSection}) = 2
number_of_parameters(::Type{TypicalSection}) = 5
inplaceness(::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{TypicalSection}) = Linear()
input_jacobian_type(::Type{TypicalSection}) = Constant()

# --- Methods --- #

function get_rates(::TypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate state rates
    return section_rhs(kh, kθ, h, θ, hdot, θdot, L, M)
end

function get_mass_matrix(::TypicalSection, q, r, p, t)
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate mass matrix
    return section_mass_matrix(m, Sθ, Iθ)
end

# --- Performance Overloads --- #

function get_state_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # return jacobian
    return section_state_jacobian(kh, kθ)
end

function get_input_jacobian(::TypicalSection)
    # return jacobian
    return section_input_jacobian()
end

# --- Unit Testing Methods --- #

function get_lhs(::TypicalSection, dq, q, r, p, t)
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # extract state rates
    dh, dθ, dhdot, dθdot = dq
    # calculate mass matrix product
    return section_lhs(m, Sθ, Iθ, dh, dθ, dhdot, dθdot)
end

# --- Internal Methods --- #

# left side of rate equations
function section_lhs(m, Sθ, Iθ, dh, dθ, dhdot, dθdot)
    SVector(dh, dθ, m*dhdot + Sθ*dθdot, Sθ*dhdot + Iθ*dθdot)
end

# right side of rate equations
function section_rhs(kh, kθ, h, θ, hdot, θdot, L, M)
    SVector(hdot, θdot, -kh*h - L, -kθ*θ + M)
end

# mass matrix (lhs jacobian wrt state rates)
function section_mass_matrix(m, Sθ, Iθ)
    @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m Sθ; 0 0 Sθ Iθ]
end

# rhs jacobian wrt states
function section_state_jacobian(kh, kθ)
    @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

# rhs jacobian wrt inputs
section_input_jacobian() = @SMatrix [0 0; 0 0; -1 0; 0 1]
