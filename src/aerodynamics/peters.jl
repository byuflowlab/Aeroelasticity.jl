"""
    Peters{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``u, \\omega,
\\dot{v}, \\dot{\\omega}`` and parameters ``a, b, a_0, \\alpha_0``
"""
struct Peters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

# --- Constructors --- #

"""
    Peters{N,TF=Float64}()

Initialize an object of type `Peters` which has `N` aerodynamic degrees of
freedom.
"""
Peters{N}() where N = Peters{N,Float64}()

function Peters{N,TF}() where {N,TF}

    b = zeros(TF, N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(big(N + n - 1))/factorial(big(N - n - 1))*
            1/factorial(big(n))^2
    end
    b[N] = (-1)^(N-1)

    c = zeros(TF, N)
    for n = 1:N
        c[n] = 2/n
    end

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[n, m] = 1/(2*n)
    end
    for m in 2:N
        n = m - 1
        D[n, m] = -1/(2*n)
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return Peters(SMatrix{N,N,TF}(A), SVector{N,TF}(b), SVector{N,TF}(c))
end

# --- Traits --- #

number_of_states(::Type{Peters{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
number_of_inputs(::Type{<:Peters}) = 4
number_of_parameters(::Type{<:Peters}) = 4

inplaceness(::Type{<:Peters}) = OutOfPlace()

rate_jacobian_type(::Type{<:Peters}) = Invariant()
state_jacobian_type(::Type{<:Peters}) = Linear()
input_jacobian_type(::Type{<:Peters}) = Nonlinear()
parameter_jacobian_type(::Type{<:Peters}) = Nonlinear()
time_gradient_type(::Type{<:Peters}) = Zeros()

# --- Methods --- #

function get_residual(model::Peters{N,TF,SV,SA}, dx, x, y, p, t) where {N,TF,SV,SA}
    # extract rates
    dλ = SVector{N}(dx)
    # extract states
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    Abar = model.A
    cbar = model.c
    # calculate rates
    return peters_residual(dλ, λ, u, ω, vdot, ωdot, a, b, Abar, cbar)
end

# --- Performance Overloads --- #

function get_rate_jacobian(model::Peters)
    # extract model constants
    Abar = model.A
    # jacobian with respect to aerodynamic states
    return peters_rate_jacobian(Abar)
end

function get_state_jacobian(model::Peters, dx, x, y, p, t)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    Abar = model.A
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return peters_state_jacobian(u, b, cbar)
end

function get_input_jacobian(model::Peters{N,TF,SV,SA}, dx, x, y, p, t) where {N,TF,SV,SA}
    # extract states
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    Abar = model.A
    cbar = model.c
    # return jacobian
    return peters_input_jacobian(λ, u, ω, a, b, Abar, cbar)
end

function get_parameter_jacobian(model::Peters{N,TF,SV,SA}, dx, x, y, p, t) where {N,TF,SV,SA}
    # extract states
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    Abar = model.A
    cbar = model.c
    # return jacobian
    return peters_parameter_jacobian(λ, u, ω, ωdot, a, b, Abar, cbar)
end

# --- Convenience Functions --- #

function set_states!(y, model::Peters; lambda)

    y .= lambda

    return y
end

function set_inputs!(y, model::Peters; u, omega, vdot, omegadot)

    y[1] = u
    y[2] = omega
    y[3] = vdot
    y[4] = omegadot

    return y
end

function set_parameters!(p, model::Peters; a, b, a0, alpha0)

    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0

    return p
end

separate_states(model::Peters, x) = (lambda = x,)

function separate_inputs(model::Peters, y)

    return (u = y[1], omega = y[2], vdot = y[3], omegadot = y[4])
end

function separate_parameters(model::Peters, p)

    return (a = p[1], b = p[2], a0 = p[3], alpha0 = p[4])
end

# --- Internal Methods for Model --- #

function peters_residual(dλ, λ, u, ω, vdot, ωdot, a, b, Abar, cbar)

    return Abar*dλ - cbar*(vdot + u*ω + (b/2-a*b)*ωdot) + u/b*λ
end

peters_rate_jacobian(Abar) = Abar

peters_state_jacobian(u, b, cbar) = u/b*Diagonal(ones(similar_type(cbar)))

function peters_input_jacobian(λ, u, ω, a, b, Abar, cbar)

    dλ_u = -cbar*ω + λ/b
    dλ_ω = -u*cbar
    dλ_vdot = -cbar
    dλ_ωdot = -(b/2-a*b)*cbar

    return hcat(dλ_u, dλ_ω, dλ_vdot, dλ_ωdot)
end

function peters_parameter_jacobian(λ, u, ω, ωdot, a, b, Abar, cbar)

    dλ_a = b*ωdot*cbar
    dλ_b = -cbar*(1/2-a)*ωdot - u/b^2*λ
    dλ_a0 = zero(cbar)
    dλ_α0 = zero(cbar)

    return hcat(dλ_a, dλ_b, dλ_a0, dλ_α0)
end

# --- Internal Methods for Couplings --- #

function peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    # circulatory load factor
    tmp1 = a0*f*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # normal force at the reference point
    N = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # axial force at the reference point
    A = -a0*ρ*b*(v + d*ω - λ0 - u*α0)^2
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*N
    # calculate loads
    return SVector(N, A, M)
end

function peters_loads_a(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3
    d = b/2 - a*b
    λ0 = 1/2 * bbar'*λ

    N = tmp1*(v + (b/2 - a*b)*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    N_a = -tmp1*b*ω - tmp2*ωdot
    A_a = 2*a0*ρ*b^2*(v + d*ω - λ0 - u*α0)*ω
    M_a = tmp2*b/2*ωdot + b*N + (b/2 + a*b)*N_a

    return SVector(N_a, A_a, M_a)
end

function peters_loads_b(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    tmp1 = a0*ρ*u*b
    tmp1_b = a0*ρ*u

    tmp2 = pi*ρ*b^3
    tmp2_b = 3*pi*ρ*b^2

    d = b/2 - a*b
    d_b = 1/2 - a

    λ0 = 1/2 * bbar'*λ

    N = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)

    N_b = tmp1_b*(v + d*ω - λ0 - u*α0) + tmp1*d_b*ω +
        tmp2_b*(vdot/b + u/b*ω - a*ωdot) + tmp2*(-vdot/b^2 - u/b^2*ω)
    A_b = -2*a0*ρ*b*(v + d*ω - λ0 - u*α0)*d_b - a0*ρ*(v + d*ω - λ0 - u*α0)^2
    M_b = -tmp2_b*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) - tmp2*(1/8 - a/2)*ωdot +
        (1/2 + a)*N + (b/2 + a*b)*N_b

    return SVector(N_b, A_b, M_b)
end

function peters_loads_ρ(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    tmp1_ρ = a0*u*b

    tmp2_ρ = pi*b^3

    d = b/2 - a*b

    λ0 = 1/2 * bbar'*λ

    N_ρ = tmp1_ρ*(v + d*ω - λ0 - u*α0) + tmp2_ρ*(vdot/b + u/b*ω - a*ωdot)
    A_ρ = -a0*b*(v + d*ω - λ0 - u*α0)^2
    M_ρ = -tmp2_ρ*(vdot/2 + u*ω + (b/8 - a*b/2)*ωdot) + (b/2 + a*b)*N_ρ

    return SVector(N_ρ, A_ρ, M_ρ)
end

function peters_loads_a0(a, b, ρ, a0, α0, bbar, u, v, ω, λ)

    tmp1_a0 = ρ*u*b

    tmp2 = pi*ρ*b^3

    d = b/2 - a*b

    λ0 = 1/2 * bbar'*λ

    N_a0 = tmp1_a0*(v + d*ω - λ0 - u*α0)
    A_a0 = -ρ*b*(v + d*ω - λ0 - u*α0)^2
    M_a0 = (b/2 + a*b)*N_a0

    return SVector(N_a0, A_a0, M_a0)
end

function peters_loads_α0(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    tmp1 = a0*ρ*u*b
    d = b/2 - a*b
    λ0 = 1/2 * bbar'*λ
    N_α0 = -tmp1*u
    A_α0 = 2*a0*ρ*b*u*(v + d*ω - λ0 - u*α0)
    M_α0 = (b/2 + a*b)*N_α0
    return SVector(N_α0, A_α0, M_α0)
end

function peters_loads_u(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    tmp1_u = a0*ρ*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # normal force at reference point
    N_u = tmp1_u*(v + d*ω - λ0 - u*α0) - tmp1*α0 + tmp2/b*ω
    # axial force at reference point
    A_u = 2*a0*ρ*b*α0*(v + d*ω - λ0 - u*α0)
    # moment at reference point
    M_u = -tmp2*ω + (b/2 + a*b)*N_u

    return SVector(N_u, A_u, M_u)
end

function peters_loads_v(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    d = b/2 - a*b
    λ0 = 1/2 * bbar'*λ
    N_v = a0*ρ*u*b
    A_v = -2*a0*ρ*b*(v + d*ω - λ0 - u*α0)
    M_v = (b/2 + a*b)*N_v
    return SVector(N_v, A_v, M_v)
end

function peters_loads_ω(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    tmp = pi*ρ*b^3
    d = b/2 - a*b
    λ0 = 1/2 * bbar'*λ
    N_ω = a0*ρ*u*b*(b/2 - a*b) + tmp*u/b
    A_ω = 2*a0*ρ*a*b^2*(v + d*ω - λ0 - u*α0)
    M_ω = -tmp*u + (b/2 + a*b)*N_ω
    return SVector(N_ω, A_ω, M_ω)
end

function peters_loads_vdot(a, b, ρ)
    tmp = pi*ρ*b^3
    N_vdot = tmp/b
    A_vdot = 0
    M_vdot = -tmp/2 + (b/2 + a*b)*N_vdot
    return SVector(N_vdot, A_vdot, M_vdot)
end

function peters_loads_ωdot(a, b, ρ)
    tmp = pi*ρ*b^3
    N_ωdot = -tmp*a
    A_ωdot = 0
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*N_ωdot
    return SVector(N_ωdot, A_ωdot, M_ωdot)
end

function peters_loads_λ(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    tmp = a0*ρ*u*b
    d = b/2 - a*b
    λ0 = 1/2 * bbar'*λ
    N_λ = -tmp/2*bbar'
    A_λ = a0*ρ*b*(v + d*ω - λ0 - u*α0)*bbar'
    M_λ = (b/2 + a*b)*N_λ
    return vcat(N_λ, A_λ, M_λ)
end

function peters_state_loads(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # normal force at reference point
    N = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*u/b*ω
    # axial force at reference point
    A = -a0*ρ*b*(v + d*ω - λ0 - u*α0)^2
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*N

    return SVector(N, A, M)
end

function peters_rate_loads(a, b, ρ, vdot, ωdot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # normal force at reference point
    N = tmp*(vdot/b - a*ωdot)
    # axial force at reference point
    A = 0
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*N

    return SVector(N, A, M)
end
