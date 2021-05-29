"""
    CoupledTypicalSection <: AbstractModel

Typical section structural model with state variables ``q = \\begin{bmatrix} h &
θ & \\dot{h} & \\dot{\\theta} \\end{bmatrix}^T``, structural parameters ``p_s =
\\begin{bmatrix} a & b & k_h & k_\\theta & m & x_\\theta & I_P \\end{bmatrix}^T``,
and aerodynamic loads ``r = \\begin{bmatrix} L & M_\\frac{1}{4} \\end{bmatrix}^T``
"""
struct CoupledTypicalSection <: AbstractModel end

# --- Traits --- #
number_of_states(::Type{CoupledTypicalSection}) = 4
number_of_inputs(::Type{CoupledTypicalSection}) = 0
number_of_parameters(::Type{CoupledTypicalSection}) = 7
inplaceness(::Type{CoupledTypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{CoupledTypicalSection}) = Constant()
state_jacobian_type(::Type{CoupledTypicalSection}) = Varying()
input_jacobian_type(::Type{CoupledTypicalSection}) = Varying()
input_dependence_type(::Type{CoupledTypicalSection}) = Linear()

# --- Methods --- #

function get_mass_matrix(::CoupledTypicalSection, q, r, p, t)
    # extract structural parameters
    U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # update mass matrix
    return section_mass_matrix(b, m, xθ, Ip)
end

function get_rates(::CoupledTypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # calculate state rates
    return section_rhs(U, ρ, a, b, kh, kθ, h, θ, hdot, θdot, L, M)
end

function get_state_jacobian(::CoupledTypicalSection, q, r, p, t)
    # extract parameters
    U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return section_state_jacobian(U, ρ, a, b, kh, kθ)
end

# TODO: Add get_parameter_jacobian

# --- Internal --- #

# left side of rate equations (used for testing)
function section_lhs(b, m, xθ, Ip, dh, dθ, dhdot, dθdot)
    SVector(dh, dθ, m*(dhdot + b*xθ*dθdot), Ip*dθdot + m*b*xθ*dhdot)
end

# right side of rate equations (used for testing)
function section_rhs(U, ρ, a, b, kh, kθ, h, θ, hdot, θdot, L, M)
    SVector(hdot, θdot, -kh*h - 2*pi*ρ*b*U^2*θ, -kθ*θ + (b/2+a*b)*2*pi*ρ*b*U^2*θ)
end

# state jacobian
function section_state_jacobian(U, ρ, a, b, kh, kθ)
    @SMatrix [0 0 1 0; 0 0 0 1; -kh -2*pi*ρ*b*U^2 0 0; 0 -kθ+(b/2+a*b)*2*pi*ρ*b*U^2 0 0]
end

# rate jacobian (mass matrix)
function section_mass_matrix(b, m, xθ, Ip)
    @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
end
