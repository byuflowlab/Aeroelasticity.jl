"""
    TypicalSection <: StructuralModel

Typical section structural model with state variables `q = [h, θ, hdot, θdot]`,
parameters `p = [a, b, kh, kθ, m, xθ, Ip]`, and aerodynamic loads `r = [L, M]`
"""
struct TypicalSection <: StructuralModel end

isinplace(::TypicalSection) = false
has_mass_matrix(::TypicalSection) = false
constant_mass_matrix(::TypicalSection) = false
linear_load_dependence(::TypicalSection) = true
defined_state_jacobian(::TypicalSection) = true
defined_load_jacobian(::TypicalSection) = true

init_mass_matrix(::TypicalSection) = zeros(4,4)

function update_mass_matrix!(::TypicalSection, M, q, r, p, t)
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # update mass matrix
    M .= section_mass_matrix(b, m, xθ, Ip)
    # return result
    return M
end

function get_rates(::TypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = u
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # extract aerodynamic loads
    L, M = r
    # calculate state rates
    return section_rhs(a, b, kh, kθ, h, θ, hdot, θdot, L, M)
end

function get_state_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return section_state_jacobian(kh, kθ)
end

function get_load_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return section_load_jacobian(a, b)
end

# TODO: Add get_parameter_jacobian

# --- Internal Functions --- #

# left side of rate equations (used for testing)
section_lhs(b, m, xθ, Ip, dh, dθ, dhdot, dθdot) = SVector(dh, dθ, m*(dhdot + b*xθ*dθdot), Ip*dθdot + m*b*xθ*dhdot)

# right side of rate equations (used for testing)
section_rhs(a, b, kh, kθ, h, θ, hdot, θdot, L, M) = SVector(hdot, θdot, -kh*h - L, -kθ*θ + M + (b/2+a*b)*L)

# state jacobian
section_state_jacobian(kh, kθ) = @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]

# load jacobian
section_load_jacobian(a, b) = @SMatrix [0 0; 0 0; -1 0; b/2+a*b 1]

# rate jacobian (mass matrix)
section_mass_matrix(b, m, xθ, Ip) = @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
