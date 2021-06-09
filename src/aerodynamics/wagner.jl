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
number_of_parameters(::Type{<:Wagner}) = 6
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
    a, b, U, ρ, a0, α0 = p
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
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    ε1 = model.eps1
    ε2 = model.eps2
    # jacobian with respect to aerodynamic states
    return wagner_state_jacobian(a, b, U, ε1, ε2)
end

function get_input_jacobian(model::Wagner, λ, d, p, t)
    # extract parameters
    a, b, U, ρ, a0, α0 = p
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
    a, b, U, ρ, a0, α0 = p
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0]
    Mrs = hcat(
        SVector(0, 0),
        SVector(0, 0),
        -wagner_loads_hddot(a, b, ρ),
        -wagner_loads_θddot(a, b, ρ)
    )
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

function get_inputs(aero::Wagner, stru::TypicalSection, u, p, t)
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = u
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # calculate aerodynamic loads
    L, M = wagner_loads(a, b, U, ρ, a0, α0, C1, C2, θ, hdot, θdot, 0, 0, λ1, λ2)
    # return portion of inputs that is not dependent on the state rates
    return SVector(θ, hdot, θdot, L, M)
end

function get_input_state_jacobian(aero::Wagner, stru::TypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # compute jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [0 1 0 0; 0 0 1 0; 0 0 0 1]
    Jra = wagner_loads_λ(a, b, U, ρ, a0)
    Jrs = hcat(
        wagner_loads_h(),
        wagner_loads_θ(a, b, U, ρ, a0, C1, C2),
        wagner_loads_hdot(a, b, U, ρ, a0, C1, C2),
        wagner_loads_θdot(a, b, U, ρ, a0, C1, C2)
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

function wagner_loads(a, b, U, ρ, a0, α0, C1, C2, θ, hdot, θdot, hddot, θddot,
    λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*U*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L = tmp1*((U*θ + hdot + d*θdot - U*α0)*ϕ0 + λ1 + λ2) + tmp2*(hddot/b + U/b*θdot - a*θddot)
    # moment at reference point
    M = -tmp2*(hddot/2 + U*θdot + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function wagner_loads_λ(a, b, U, ρ, a0)
    tmp1 = a0*ρ*U*b
    tmp2 = (b/2 + a*b)*tmp1
    return @SMatrix [tmp1 tmp1; tmp2 tmp2]
end
wagner_loads_λdot() = @SMatrix [0 0; 0 0]
wagner_loads_h() = SVector(0, 0)
function wagner_loads_θ(a, b, U, ρ, a0, C1, C2)
    ϕ0 = 1 - C1 - C2
    L_θ = a0*ρ*b*U^2*ϕ0
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end
function wagner_loads_hdot(a, b, U, ρ, a0, C1, C2)
    ϕ0 = 1 - C1 - C2
    L_hdot = a0*ρ*U*b*ϕ0
    M_hdot = (b/2 + a*b)*L_hdot
    return SVector(L_hdot, M_hdot)
end
function wagner_loads_θdot(a, b, U, ρ, a0, C1, C2)
    tmp = pi*ρ*b^3
    ϕ0 = 1 - C1 - C2
    L_θdot = a0*ρ*U*b*(b/2 - a*b)*ϕ0 + tmp*U/b
    M_θdot = -tmp*U + (b/2 + a*b)*L_θdot
    return SVector(L_θdot, M_θdot)
end
function wagner_loads_hddot(a, b, ρ)
    tmp1 = pi*ρ*b^3
    tmp2 = b/2 + a*b
    L_hddot = tmp1/b
    M_hddot = -tmp1/2 + tmp2*L_hddot
    return SVector(L_hddot, M_hddot)
end
function wagner_loads_θddot(a, b, ρ)
    tmp1 = pi*ρ*b^3
    tmp2 = b/2 + a*b
    L_θddot = -tmp1*a
    M_θddot = -tmp1*(b/8 - a*b/2) + tmp2*L_θddot
    return SVector(L_θddot, M_θddot)
end
