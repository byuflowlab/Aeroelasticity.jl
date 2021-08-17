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
mass_matrix_type(::Type{<:Peters}) = Constant()
state_jacobian_type(::Type{<:Peters}) = Linear()
input_jacobian_type(::Type{<:Peters}) = Nonlinear()

# --- Methods --- #

function get_rates(model::Peters{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract inputs
    u, ω, vdot, ωdot = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return peters_rhs(a, b, cbar, u, ω, vdot, ωdot, λ)
end

get_mass_matrix(model::Peters) = model.A

# --- Performance Overloads --- #

function get_state_jacobian(model::Peters, λ, d, p, t)
    # extract inputs
    u, ω, vdot, ωdot = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return peters_state_jacobian(b, cbar, u)
end

function get_input_jacobian(model::Peters{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract inputs
    u, ω, vdot, ωdot = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return peters_input_jacobian(a, b, cbar, u, ω, λ)
end

# --- Unit Testing Methods --- #

get_lhs(model::Peters, dλ, λ, d, p, t) = model.A*dλ

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


# --- Internal Methods --- #

peters_lhs(Abar, dλ) = Abar*dλ

peters_rhs(a, b, cbar, u, ω, vdot, ωdot, λ) = cbar*(vdot + u*ω + (b/2-a*b)*ωdot) - u/b*λ

peters_mass_matrix(Abar) = Abar

peters_state_jacobian(b, cbar, u) = -u/b*Diagonal(one.(cbar))

function peters_input_jacobian(a, b, cbar, u, ω, λ)
    return hcat(cbar*ω - λ/b, u*cbar, cbar, (b/2-a*b)*cbar)
end

function peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
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
    # lift at reference point
    L = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L

    return SVector(L, M)
end

function peters_rate_loads(a, b, ρ, vdot, ωdot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(vdot/b - a*ωdot)
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function peters_loads_λ(a, b, ρ, a0, bbar, u)
    tmp = a0*ρ*u*b
    L_λ = -tmp/2*bbar'
    M_λ = (b/2 + a*b)*L_λ
    return vcat(L_λ, M_λ)
end

function peters_loads_λdot(cbar)
    tmp = zero(cbar)'
    return vcat(tmp, tmp)
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
    # lift at reference point
    L_u = tmp1_u*(v + d*ω - λ0 - u*α0) - tmp1*α0 + tmp2/b*ω
    # moment at reference point
    M_u = -tmp2*ω + (b/2 + a*b)*L_u

    return SVector(L_u, M_u)
end

function peters_loads_v(a, b, ρ, a0, u)
    # lift at reference point
    L_v = a0*ρ*u*b
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

function peters_loads_ω(a, b, ρ, a0, u)
    tmp = pi*ρ*b^3
    L_ω = a0*ρ*u*b*(b/2 - a*b) + tmp*u/b
    M_ω = -tmp*u + (b/2 + a*b)*L_ω
    return SVector(L_ω, M_ω)
end

peters_loads_udot() = SVector(0, 0)

function peters_loads_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_vdot = tmp/b
    # moment at reference point
    M_vdot = -tmp/2 + (b/2 + a*b)*L_vdot

    return SVector(L_vdot, M_vdot)
end

function peters_loads_ωdot(a, b, ρ)
    tmp = pi*ρ*b^3
    L_ωdot = -tmp*a
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_ωdot
    return SVector(L_ωdot, M_ωdot)
end

function peters_loads_θ(a, b, ρ, a0, u)
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end
