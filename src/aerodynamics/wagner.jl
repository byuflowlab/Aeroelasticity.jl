"""
    Wagner{TF} <: AbstractModel

Aerodynamic model based on Wagner's function with state variables ``\\lambda_1,
\\lambda_2``, inputs ``u, v, \\omega``, and parameters ``a, b, a_0, \\alpha_0``
"""
struct Wagner{TF} <: AbstractModel
    C1::TF
    C2::TF
    eps1::TF
    eps2::TF
end

# --- Constructors --- #

"""
    Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

Initialize an object of type [`Wagner`](@ref)
"""
function Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)
    return Wagner(C1, C2, eps1, eps2)
end

# --- Traits --- #

number_of_states(::Type{<:Wagner}) = 2
number_of_inputs(::Type{<:Wagner}) = 3
number_of_parameters(::Type{<:Wagner}) = 4
inplaceness(::Type{<:Wagner}) = OutOfPlace()
mass_matrix_type(::Type{<:Wagner}) = Identity()
state_jacobian_type(::Type{<:Wagner}) = Linear()
input_jacobian_type(::Type{<:Wagner}) = Nonlinear()

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

# --- Internal Methods --- #

function wagner_rates(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)
    λ1dot = -ε1*u/b*λ1 + C1*ε1*u/b*(v + (b/2-a*b)*ω - u*α0)
    λ2dot = -ε2*u/b*λ2 + C2*ε2*u/b*(v + (b/2-a*b)*ω - u*α0)
    return SVector(λ1dot, λ2dot)
end

wagner_state_jacobian(a, b, ε1, ε2, u) = @SMatrix [-ε1*u/b 0; 0 -ε2*u/b]

function wagner_input_jacobian(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)

    tmp1 = v/b + (1/2-a)*ω - 2*α0*u/b
    λ1dot_u = -ε1/b*λ1 + C1*ε1*tmp1
    λ2dot_u = -ε2/b*λ2 + C2*ε2*tmp1

    tmp2 = u/b
    λ1dot_v = C1*ε1*tmp2
    λ2dot_v = C2*ε2*tmp2

    tmp3 = u*(1/2-a)
    λ1dot_ω = C1*ε1*tmp3
    λ2dot_ω = C2*ε2*tmp3

    return @SMatrix [λ1dot_u λ1dot_v λ1dot_ω; λ2dot_u λ2dot_v λ2dot_ω]
end

function wagner_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function wagner_state_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L

    return SVector(L, M)
end

function wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(vdot/b - a*ωdot)
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function wagner_loads_λ(a, b, ρ, a0, u)
    tmp1 = a0*ρ*u*b
    tmp2 = (b/2 + a*b)*tmp1
    return @SMatrix [tmp1 tmp1; tmp2 tmp2]
end
wagner_loads_λdot() = @SMatrix [0 0; 0 0]

function wagner_loads_u(a, b, ρ, a0, α0, C1, C2, u, v, ωdot, λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    tmp1_u = a0*ρ*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_u = tmp1_u*((v + d*ωdot - u*α0)*ϕ0 + λ1 + λ2) - tmp1*α0*ϕ0 + tmp2/b*ωdot
    # moment at reference point
    M_u = -tmp2*ωdot + (b/2 + a*b)*L_u

    return SVector(L_u, M_u)
end

function wagner_loads_v(a, b, ρ, a0, C1, C2, u)
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_v = a0*ρ*u*b*ϕ0
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

function wagner_loads_ω(a, b, ρ, a0, C1, C2, u)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_ω = tmp1*d*ϕ0 + tmp2*u/b
    # moment at reference point
    M_ω = -tmp2*u + (b/2 + a*b)*L_ω

    return SVector(L_ω, M_ω)
end

wagner_loads_udot() = SVector(0, 0)

function wagner_loads_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_vdot = tmp/b
    # moment at reference point
    M_vdot = -tmp/2 + (b/2 + a*b)*L_vdot

    return SVector(L_vdot, M_vdot)
end

function wagner_loads_ωdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_ωdot = -a*tmp
    # moment at reference point
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_ωdot

    return SVector(L_ωdot, M_ωdot)
end

wagner_loads_h() = SVector(0, 0)

function wagner_loads_θ(a, b, ρ, a0, C1, C2, u)
    ϕ0 = 1 - C1 - C2
    L_θ = a0*ρ*b*u^2*ϕ0
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end
