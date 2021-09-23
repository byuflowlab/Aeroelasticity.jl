"""
    DeformableAirfoil{N,TF,SV,SA} <: AbstractModel

Deformable airfoil model with `N` state variables, inputs ``u, \\omega,
\\dot{v}, \\dot{\\omega}`` and parameters ``a, b, a_0, \\alpha_0``
"""
struct DeformableAirfoil{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

# --- Constructors --- #

"""
    DeformableAirfoil{N,TF=Float64}()

Initialize an object of type `DeformableAirfoil` which has `N` aerodynamic degrees of
freedom.
"""
DeformableAirfoil{N}() where N = DeformableAirfoil{N,Float64}()

function DeformableAirfoil{N,TF}() where {N,TF}

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

    return DeformableAirfoil(SMatrix{N,N,TF}(A), SVector{N,TF}(b), SVector{N,TF}(c))
end

# --- Traits --- #

number_of_states(::Type{DeformableAirfoil{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
number_of_inputs(::Type{<:DeformableAirfoil}) = 4
number_of_parameters(::Type{<:DeformableAirfoil}) = 4

inplaceness(::Type{<:DeformableAirfoil}) = OutOfPlace()

rate_jacobian_type(::Type{<:DeformableAirfoil}) = Invariant()
state_jacobian_type(::Type{<:DeformableAirfoil}) = Linear()
input_jacobian_type(::Type{<:DeformableAirfoil}) = Nonlinear()
parameter_jacobian_type(::Type{<:DeformableAirfoil}) = Nonlinear()
time_gradient_type(::Type{<:DeformableAirfoil}) = Zeros()

# --- Methods --- #

function get_residual(model::DeformableAirfoil{N,TF,SV,SA}, dx, x, y, p, t) where {N,TF,SV,SA}
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

# TODO

# --- Convenience Functions --- #

# TODO

# --- Internal Methods for Model --- #

# Come from Peters Finite State Model

# --- Internal Methods for Couplings --- #

function generalized_motion(u, v, ω)

    u0 = u
    w0 = v + u0*θ/b
    w1 = ω

end

function deformable_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    L = 2*pi*ρ*b*f*u*(v - λ0 - u*α0) - pi*ρ*b*u*ω - pi*ρ*b^2*vdot

    2*pi*ρ*b^2*u*ω

    L = a0*ρ*u*b*(v + (b/2 - a*b)*ω - λ0 - u*α0) + pi*ρ*b^2*(vdot) + pi*ρ*b*udot*h0*h1 + 2*pi*ρ*u^2*sum(1:N .* hn)
    L = a0*ρ*u*b*(v + (b/2 - a*b)*ω - λ0 - u*α0) + pi*ρ*b^3*(vdot/b + u/b*ω - a*ωdot)



    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    pi*ρ*b^2*()
    L = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L
    # calculate loads
    return SVector(L, M)
end

function deformable_loads_a(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3

    λ0 = 1/2 * bbar'*λ

    L = tmp1*(v + (b/2 - a*b)*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    L_a = -tmp1*b*ω - tmp2*ωdot

    M_a = tmp2*b/2*ωdot + b*L + (b/2 + a*b)*L_a

    return SVector(L_a, M_a)
end

function deformable_loads_b(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    tmp1 = a0*ρ*u*b
    tmp1_b = a0*ρ*u

    tmp2 = pi*ρ*b^3
    tmp2_b = 3*pi*ρ*b^2

    d = b/2 - a*b
    d_b = 1/2 - a

    λ0 = 1/2 * bbar'*λ

    L = tmp1*(v + d*ω - λ0 - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)

    L_b = tmp1_b*(v + d*ω - λ0 - u*α0) + tmp1*d_b*ω +
        tmp2_b*(vdot/b + u/b*ω - a*ωdot) + tmp2*(-vdot/b^2 - u/b^2*ω)
    M_b = -tmp2_b*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) - tmp2*(1/8 - a/2)*ωdot +
        (1/2 + a)*L + (b/2 + a*b)*L_b

    return SVector(L_b, M_b)
end

function deformable_loads_ρ(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    tmp1_ρ = a0*u*b

    tmp2_ρ = pi*b^3

    d = b/2 - a*b

    λ0 = 1/2 * bbar'*λ

    L_ρ = tmp1_ρ*(v + d*ω - λ0 - u*α0) + tmp2_ρ*(vdot/b + u/b*ω - a*ωdot)

    M_ρ = -tmp2_ρ*(vdot/2 + u*ω + (b/8 - a*b/2)*ωdot) + (b/2 + a*b)*L_ρ

    return SVector(L_ρ, M_ρ)
end

function deformable_loads_a0(a, b, ρ, a0, α0, bbar, u, v, ω, λ)

    tmp1_a0 = ρ*u*b

    tmp2 = pi*ρ*b^3

    d = b/2 - a*b

    λ0 = 1/2 * bbar'*λ

    L_a0 = tmp1_a0*(v + d*ω - λ0 - u*α0)
    M_a0 = (b/2 + a*b)*L_a0

    return SVector(L_a0, M_a0)
end

function deformable_loads_α0(a, b, ρ, a0, u)

    tmp1 = a0*ρ*u*b

    L_α0 = -tmp1*u

    M_α0 = (b/2 + a*b)*L_α0

    return SVector(L_α0, M_α0)
end

function deformable_loads_u(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
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

function deformable_loads_v(a, b, ρ, a0, u)
    # lift at reference point
    L_v = a0*ρ*u*b
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

function deformable_loads_ω(a, b, ρ, a0, u)
    tmp = pi*ρ*b^3
    L_ω = a0*ρ*u*b*(b/2 - a*b) + tmp*u/b
    M_ω = -tmp*u + (b/2 + a*b)*L_ω
    return SVector(L_ω, M_ω)
end

function deformable_loads_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_vdot = tmp/b
    # moment at reference point
    M_vdot = -tmp/2 + (b/2 + a*b)*L_vdot

    return SVector(L_vdot, M_vdot)
end

function deformable_loads_ωdot(a, b, ρ)
    tmp = pi*ρ*b^3
    L_ωdot = -tmp*a
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_ωdot
    return SVector(L_ωdot, M_ωdot)
end

function deformable_loads_λ(a, b, ρ, a0, bbar, u)
    tmp = a0*ρ*u*b
    L_λ = -tmp/2*bbar'
    M_λ = (b/2 + a*b)*L_λ
    return vcat(L_λ, M_λ)
end

function deformable_state_loads(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
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

function deformable_rate_loads(a, b, ρ, vdot, ωdot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(vdot/b - a*ωdot)
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end
