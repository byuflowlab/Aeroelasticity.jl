"""
    TypicalSection <: AbstractModel

Typical section structural model with state variables ``q = \\begin{bmatrix} h &
θ & \\dot{h} & \\dot{\\theta} \\end{bmatrix}^T``, structural parameters ``p_s =
\\begin{bmatrix} a & b & k_h & k_\\theta & m & x_\\theta & I_P \\end{bmatrix}^T``,
and aerodynamic loads ``r = \\begin{bmatrix} L & M_\\frac{1}{4} \\end{bmatrix}^T``
"""
struct TypicalSection <: AbstractModel end

# --- Traits --- #
number_of_states(::Type{TypicalSection}) = 4
number_of_inputs(::Type{TypicalSection}) = 2
number_of_parameters(::Type{TypicalSection}) = 7
inplaceness(::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{TypicalSection}) = Constant()
state_jacobian_type(::Type{TypicalSection}) = Varying()
input_jacobian_type(::Type{TypicalSection}) = Varying()
input_dependence_type(::Type{TypicalSection}) = Linear()

# --- Methods --- #

function get_mass_matrix(::TypicalSection, q, r, p, t)
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # update mass matrix
    return section_mass_matrix(b, m, xθ, Ip)
end

function get_rates(::TypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # calculate state rates
    return section_rhs(a, b, kh, kθ, h, θ, hdot, θdot, L, M)
end

function get_state_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return section_state_jacobian(kh, kθ)
end

function get_input_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return section_load_jacobian(a, b)
end

# TODO: Add get_parameter_jacobian

# --- Internal --- #

# left side of rate equations (used for testing)
function section_lhs(b, m, xθ, Ip, dh, dθ, dhdot, dθdot)
    SVector(dh, dθ, m*(dhdot + b*xθ*dθdot), Ip*dθdot + m*b*xθ*dhdot)
end

# right side of rate equations (used for testing)
function section_rhs(a, b, kh, kθ, h, θ, hdot, θdot, L, M)
    SVector(hdot, θdot, -kh*h - L, -kθ*θ + M + (b/2+a*b)*L)
end

# state jacobian
function section_state_jacobian(kh, kθ)
    @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

# load jacobian
section_load_jacobian(a, b) = @SMatrix [0 0; 0 0; -1 0; b/2+a*b 1]

# rate jacobian (mass matrix)
function section_mass_matrix(b, m, xθ, Ip)
    @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
end
