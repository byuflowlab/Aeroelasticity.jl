"""
    Wagner{TF} <: AbstractModel

Aerodynamic model based on Wagner's function with state variables ``d =
\\begin{bmatrix} \\lambda_1 & \\lambda_2 \\end{bmatrix}^T``, inputs ``d =
\\begin{bmatrix} \\theta & \\dot{h} & \\dot{\\theta} \\end{bmatrix}^T``
and parameters ``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``
"""
struct Wagner{TF} <: AbstractModel
    C1::TF
    C2::TF
    eps1::TF
    eps2::TF
end

# --- Constructors --- #
function Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)
    return Wagner(C1, C2, eps1, eps2)
end

# --- Traits --- #
number_of_states(::Type{<:Wagner}) = 2
number_of_inputs(::Type{<:Wagner}) = 3
number_of_parameters(::Type{<:Wagner}) = 4
inplaceness(::Type{<:Wagner}) = OutOfPlace()
mass_matrix_type(::Type{<:Wagner}) = Identity()
state_jacobian_type(::Type{<:Wagner}) = Varying()
input_jacobian_type(::Type{<:Wagner}) = Varying()
input_dependence_type(::Type{<:Wagner}) = Linear()

# --- Methods --- #

function get_rates(model::Wagner, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states
    λ1, λ2 = λ
    # extract structural deflections
    θ, hdot, θdot = d
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # calculate rates
    return wagner_rhs(a, b, U, C1, C2, ε1, ε2, θ, hdot, θdot, λ1, λ2)
end

function get_state_jacobian(model::Wagner, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    ε1 = model.eps1
    ε2 = model.eps2
    # jacobian with respect to aerodynamic states
    return wagner_state_jacobian(a, b, U, ε1, ε2)
end

function get_input_jacobian(model::Wagner, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # return jacobian
    return wagner_input_jacobian(a, b, U, C1, C2, ε1, ε2)
end

# TODO: Add parameter jacobian

# --- Coupled Model Properties --- #

# traits
inplaceness(::Type{<:Wagner}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Varying()
state_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Varying()

function get_input_mass_matrix(aero::Wagner, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0]
    Mrs = hcat(
        SVector(0, 0),
        SVector(0, 0),
        -wagner_loads_hddot(b, ρ),
        -wagner_loads_θddot(a, b, ρ)
    )
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

function get_inputs(aero::Wagner, stru::TypicalSection, u, p, t)
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = u
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # calculate aerodynamic loads
    L, M = wagner_loads(a, b, U, ρ, C1, C2, θ, hdot, θdot, 0, 0, λ1, λ2)
    # return portion of inputs that is not dependent on the state rates
    return SVector(θ, hdot, θdot, L, M)
end

function get_input_state_jacobian(aero::Wagner, stru::TypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # compute jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [0 1 0 0; 0 0 1 0; 0 0 0 1]
    Jra = wagner_loads_λ(b, U, ρ)
    Jrs = hcat(
        wagner_loads_h(),
        wagner_loads_θ(b, U, ρ, C1, C2),
        wagner_loads_hdot(b, U, ρ, C1, C2),
        wagner_loads_θdot(a, b, U, ρ, C1, C2)
        )
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# TODO: Parameter jacobian

# --- Internal Methods --- #

function wagner_rates(a, b, U, C1, C2, ε1, ε2, θ, hdot, θdot, λ1, λ2)
    λ1dot = -ε1*U/b*λ1 + C1*ε1*U/b*(hdot + (1/2-a)*b*θdot + U*θ)
    λ2dot = -ε2*U/b*λ2 + C2*ε2*U/b*(hdot + (1/2-a)*b*θdot + U*θ)
    return SVector(λ1dot, λ2dot)
end

wagner_state_jacobian(a, b, U, ε1, ε2) = @SMatrix [-ε1*U/b 0; 0 -ε2*U/b]

function wagner_input_jacobian(a, b, U, C1, C2, ε1, ε2)
    tmp1 = C1*ε1*U
    tmp2 = C2*ε2*U
    return @SMatrix [tmp1*U/b tmp1/b tmp1/2-a*tmp1; tmp2*U/b tmp2/b tmp2/2-a*tmp2]
end

function wagner_loads(a, b, U, ρ, C1, C2, θ, hdot, θdot, hddot, θddot, λ1, λ2)
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # downwash at the 3/4 chord
    wt = hdot + (1/2 - a)*b*θdot + U*θ
    # lift
    L = 2*pi*ρ*b*U*(wt*ϕ0 + λ1 + λ2) + pi*ρ*b^2*(hddot + U*θdot - a*b*θddot)
    # moment at the 1/4 chord
    M = -pi*ρ*b^3*(hddot/2 + U*θdot + b*(1/8 - a/2)*θddot)
    # return load
    return SVector(L, M)
end

wagner_loads_λ(b, U, ρ) = 2*pi*ρ*b*U*(@SMatrix [1 1; 0 0])
wagner_loads_λdot() = @SMatrix [0 0; 0 0]
wagner_loads_h() = SVector(0, 0)
function wagner_loads_θ(b, U, ρ, C1, C2)
    ϕ0 = 1 - C1 - C2
    return SVector(2*pi*ρ*b*U^2*ϕ0, 0)
end
function wagner_loads_hdot(b, U, ρ, C1, C2)
    ϕ0 = 1 - C1 - C2
    return SVector(2*pi*ρ*b*U*ϕ0, 0)
end
function wagner_loads_θdot(a, b, U, ρ, C1, C2)
    ϕ0 = 1 - C1 - C2
    return SVector(pi*ρ*b^2*U*(ϕ0 - 2*ϕ0*a + 1), -pi*ρ*b^3*U)
end
wagner_loads_hddot(b, ρ) = SVector(pi*ρ*b^2, -pi/2*ρ*b^3)
wagner_loads_θddot(a, b, ρ) = SVector(-pi*ρ*a*b^3, -pi/8*ρ*b^4*(1 - 4*a))
