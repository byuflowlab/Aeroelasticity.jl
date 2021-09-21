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

rate_jacobian_type(::Type{<:Wagner}) = Identity()
state_jacobian_type(::Type{<:Wagner}) = Linear()
input_jacobian_type(::Type{<:Wagner}) = Nonlinear()
parameter_jacobian_type(::Type{<:Wagner}) = Nonlinear()
time_gradient_type(::Type{<:Wagner}) = Zeros()

# --- Methods --- #

function get_residual(model::Wagner, dx, x, y, p, t)
    # extract rates
    dλ1, dλ2 = dx
    # extract states
    λ1, λ2 = x
    # extract inputs
    u, v, ω = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # calculate rates
    return wagner_residual(dλ1, dλ2, λ1, λ2, u, v, ω, a, b, α0, C1, C2, ε1, ε2)
end

# --- Performance Overloads --- #

function get_state_jacobian(model::Wagner, dx, x, y, p, t)
    # extract inputs
    u, v, ω = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    ε1 = model.eps1
    ε2 = model.eps2
    # jacobian with respect to aerodynamic states
    return wagner_state_jacobian(u, a, b, ε1, ε2)
end

function get_input_jacobian(model::Wagner, dx, x, y, p, t)
    # extract states
    λ1, λ2 = x
    # extract inputs
    u, v, ω = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # return jacobian
    return wagner_input_jacobian(λ1, λ2, u, v, ω, a, b, α0, C1, C2, ε1, ε2)
end

function get_parameter_jacobian(model::Wagner, dx, x, y, p, t)
    # extract states
    λ1, λ2 = x
    # extract inputs
    u, v, ω = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # return jacobian
    return wagner_parameter_jacobian(λ1, λ2, u, v, ω, b, α0, C1, C2, ε1, ε2)
end

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

function wagner_residual(dλ1, dλ2, λ1, λ2, u, v, ω, a, b, α0, C1, C2, ε1, ε2)
    r1 = dλ1 + ε1*u/b*λ1 - C1*ε1*u/b*(v + (b/2-a*b)*ω - u*α0)
    r2 = dλ2 + ε2*u/b*λ2 - C2*ε2*u/b*(v + (b/2-a*b)*ω - u*α0)
    return SVector(r1, r2)
end

wagner_state_jacobian(u, a, b, ε1, ε2) = @SMatrix [ε1*u/b 0; 0 ε2*u/b]

function wagner_input_jacobian(λ1, λ2, u, v, ω, a, b, α0, C1, C2, ε1, ε2)

    tmp1 = v/b + (1/2-a)*ω - 2*α0*u/b
    λ1dot_u = ε1/b*λ1 - C1*ε1*tmp1
    λ2dot_u = ε2/b*λ2 - C2*ε2*tmp1

    tmp2 = u/b
    λ1dot_v = -C1*ε1*tmp2
    λ2dot_v = -C2*ε2*tmp2

    tmp3 = u*(1/2-a)
    λ1dot_ω = -C1*ε1*tmp3
    λ2dot_ω = -C2*ε2*tmp3

    return @SMatrix [λ1dot_u λ1dot_v λ1dot_ω; λ2dot_u λ2dot_v λ2dot_ω]
end

function wagner_parameter_jacobian(λ1, λ2, u, v, ω, b, α0, C1, C2, ε1, ε2)

    λ1dot_a = C1*ε1*u*ω
    λ1dot_b = -ε1*u/b^2*λ1 + C1*ε1/b^2*(u*v - u^2*α0)
    λ1dot_a0 = 0
    λ1dot_α0 = C1*ε1*u^2/b

    λ2dot_a = C2*ε2*u*ω
    λ2dot_b = -ε2*u/b^2*λ2 + C2*ε2/b^2*(u*v - u^2*α0)
    λ2dot_a0 = 0
    λ2dot_α0 = C2*ε2*u^2/b

    return @SMatrix [
        λ1dot_a λ1dot_b λ1dot_a0 λ1dot_α0;
        λ2dot_a λ2dot_b λ2dot_a0 λ2dot_α0;
        ]
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
    # normal force at reference point
    N = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # axial force at reference point
    A = -a0*ρ*b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)^2
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*N

    return SVector(N, A, M)
end

function wagner_loads_λ(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    d = b/2 - a*b
    ϕ0 = 1 - C1 - C2
    N_λ = a0*ρ*u*b
    A_λ = -2*a0*ρ*b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)
    M_λ = (b/2 + a*b)*N_λ
    return @SMatrix [N_λ N_λ; A_λ A_λ; M_λ M_λ]
end
wagner_loads_λdot() = @SMatrix [0 0; 0 0; 0 0]

wagner_loads_udot() = SVector(0, 0, 0)

function wagner_loads_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # normal force at reference point
    N_vdot = tmp/b
    # axial force at reference point
    A_vdot = 0
    # moment at reference point
    M_vdot = -tmp/2 + (b/2 + a*b)*N_vdot

    return SVector(N_vdot, A_vdot, M_vdot)
end

function wagner_loads_ωdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # normal force at reference point
    N_ωdot = -a*tmp
    # axial force at reference point
    A_ωdot = 0.0
    # moment at reference point
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*N_ωdot

    return SVector(N_ωdot, A_ωdot, M_ωdot)
end

function wagner_loads_a(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)

    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3

    ϕ0 = 1 - C1 - C2

    N = tmp1*((v + (b/2 - a*b)*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    N_a = -tmp1*b*ω*ϕ0 - tmp2*ωdot

    A_a = 2*a0*ρ*b^2*ω*ϕ0*((v + (b/2 - a*b)*ω - u*α0)*ϕ0 + λ1 + λ2)

    M_a = tmp2*b/2*ωdot + b*N + (b/2 + a*b)*N_a

    return SVector(N_a, A_a, M_a)
end

function wagner_loads_b(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)

    tmp1 = a0*ρ*u*b
    tmp1_b = a0*ρ*u

    tmp2 = pi*ρ*b^3
    tmp2_b = 3*pi*ρ*b^2

    d = b/2 - a*b
    d_b = 1/2 - a

    ϕ0 = 1 - C1 - C2

    N = tmp1*((v + (b/2 - a*b)*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*(vdot/b + u/b*ω - a*ωdot)

    N_b = tmp1_b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp1*d_b*ω*ϕ0 +
        tmp2_b*(vdot/b + u/b*ω - a*ωdot) + tmp2*(-vdot/b^2 - u/b^2*ω)
    A_b = -a0*ρ*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)^2 - 2*a0*ρ*b*d_b*ω*ϕ0*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)
    M_b = -tmp2_b*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) - tmp2*(1/8 - a/2)*ωdot +
        (1/2 + a)*N + (b/2 + a*b)*N_b

    return SVector(N_b, A_b, M_b)
end

function wagner_loads_ρ(a, b, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)

    tmp1_ρ = a0*u*b

    tmp2_ρ = pi*b^3

    d = b/2 - a*b

    ϕ0 = 1 - C1 - C2

    N_ρ = tmp1_ρ*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2_ρ*(vdot/b + u/b*ω - a*ωdot)
    A_ρ = -a0*b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)^2
    M_ρ = -tmp2_ρ*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*N_ρ

    return SVector(N_ρ, A_ρ, M_ρ)
end

function wagner_loads_a0(a, b, ρ, α0, C1, C2, u, v, ω, λ1, λ2)

    tmp1_a0 = ρ*u*b

    tmp2 = pi*ρ*b^3

    d = b/2 - a*b

    ϕ0 = 1 - C1 - C2

    N_a0 = tmp1_a0*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)
    A_a0 = -ρ*b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)^2
    M_a0 = (b/2 + a*b)*N_a0

    return SVector(N_a0, A_a0, M_a0)
end

function wagner_loads_α0(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)

    tmp1 = a0*ρ*u*b

    d = b/2 - a*b

    ϕ0 = 1 - C1 - C2

    N_α0 = -tmp1*u*ϕ0
    A_α0 = 2*a0*ρ*b*u*ϕ0*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)
    M_α0 = (b/2 + a*b)*N_α0

    return SVector(N_α0, A_α0, M_α0)
end

function wagner_loads_u(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    tmp1 = a0*ρ*u*b
    tmp1_u = a0*ρ*b
    tmp2 = pi*ρ*b^3
    d = b/2 - a*b
    ϕ0 = 1 - C1 - C2
    N_u = tmp1_u*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) - tmp1*α0*ϕ0 + tmp2/b*ω
    A_u = 2*a0*ρ*b*α0*ϕ0*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)
    M_u = -tmp2*ω + (b/2 + a*b)*N_u
    return SVector(N_u, A_u, M_u)
end

function wagner_loads_v(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    ϕ0 = 1 - C1 - C2
    d = b/2 - a*b
    N_v = a0*ρ*u*b*ϕ0
    A_v = -2*a0*ρ*b*ϕ0*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)
    M_v = (b/2 + a*b)*N_v

    return SVector(N_v, A_v, M_v)
end

function wagner_loads_ω(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3
    d = b/2 - a*b
    ϕ0 = 1 - C1 - C2
    N_ω = tmp1*d*ϕ0 + tmp2*u/b
    A_ω = -2*a0*ρ*b*d*ϕ0*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)
    M_ω = -tmp2*u + (b/2 + a*b)*N_ω

    return SVector(N_ω, A_ω, M_ω)
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
    # normal force at reference point
    N = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*u/b*ω
    # axial force at reference point
    A = -a0*ρ*b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)^2
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*N

    return SVector(N, A, M)
end

function wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # normal force at reference point
    N = tmp*(vdot/b - a*ωdot)
    # axial force at reference point
    A = -a0*ρ*b*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2)^2
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*N

    return SVector(N, A, M)
end
