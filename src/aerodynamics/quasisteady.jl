"""
    QuasiSteady{Order} <: NoStateModel

2D quasi-steady aerodynamic model with parameters ``a, b, a_0, \\alpha_0``.
"""
struct QuasiSteady{Order} <: NoStateModel end

"""
    Steady <: NoStateModel

2D steady aerodynamic model with parameters ``a, b, a_0, \\alpha_0``.
"""
const Steady = QuasiSteady{0}

# --- Constructors --- #

"""
    Steady()

Initialize an object of type [`Steady`](@ref) which represents a 2D
steady aerodynamic model.
"""
Steady()

"""
    QuasiSteady()

Initialize an object of type [`QuasiSteady`](@ref) which represents a 2D
quasi-steady aerodynamic model.
"""
QuasiSteady() = QuasiSteady{2}()

# --- Traits --- #

number_of_parameters(::Type{<:QuasiSteady}) = 4

# --- Convenience Methods --- #

function set_parameters!(p, model::QuasiSteady; a, b, a0, alpha0)

    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0

    return p
end

function separate_parameters(model::QuasiSteady, p)

    return (a=p[1], b=p[2], a0=p[3], alpha0=p[4])
end

# --- Internal Methods --- #

# steady state loads
function quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # lift at reference point
    L = a0*ρ*b*u*(v - u*α0)
    # moment at reference point
    M = (b/2 + a*b)*L

    return SVector(L, M)
end

function quasisteady0_loads_a(a, b, ρ, a0, α0, u, v)

    L = a0*ρ*b*u*(v - u*α0)

    L_a = 0

    M_a = b*L

    return SVector(L_a, M_a)
end

function quasisteady0_loads_b(a, b, ρ, a0, α0, u, v)

    L = a0*ρ*b*u*(v - u*α0)

    L_b = a0*ρ*u*(v - u*α0)

    M_b = (1/2 + a)*L + (b/2 + a*b)*L_b

    return SVector(L_b, M_b)
end

function quasisteady0_loads_ρ(a, b, ρ, a0, α0, u, v)
    # lift at reference point
    L_ρ = a0*b*u*(v - u*α0)
    # moment at reference point
    M_ρ = (b/2 + a*b)*L_ρ

    return SVector(L_ρ, M_ρ)
end

function quasisteady0_loads_a0(a, b, ρ, α0, u, v)
    # lift at reference point
    L_a0 = ρ*b*u*(v - u*α0)
    # moment at reference point
    M_a0 = (b/2 + a*b)*L_a0

    return SVector(L_a0, M_a0)
end

function quasisteady0_loads_α0(a, b, ρ, a0, α0, u, v)
    # lift at reference point
    L_α0 = -a0*ρ*b*u^2
    # moment at reference point
    M_α0 = (b/2 + a*b)*L_α0

    return SVector(L_α0, M_α0)
end

function quasisteady0_u(a, b, ρ, a0, v)
    L_u = a0*ρ*b*v
    M_u = (b/2 + a*b)*L_u
    return SVector(L_u, M_u)
end

function quasisteady0_v(a, b, ρ, a0, u)
    L_v = a0*ρ*b*u
    M_v = (b/2 + a*b)*L_v
    return SVector(L_v, M_v)
end

# circulatory loads + added mass effects, neglecting accelerations
function quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L = tmp1*(v + d*ω - u*α0) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L

    return SVector(L, M)
end

function quasisteady1_loads_a(a, b, ρ, a0, α0, u, v, ω)

    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3

    L = tmp1*(v + (b/2 - a*b)*ω - u*α0) + tmp2*u/b*ω
    L_a = -tmp1*b*ω

    M_a = b*L + (b/2 + a*b)*L_a

    return SVector(L_a, M_a)
end


function quasisteady1_loads_b(a, b, ρ, a0, α0, u, v, ω)

    tmp1 = a0*ρ*u*b
    tmp1_b = a0*ρ*u

    tmp2 = pi*ρ*b^3
    tmp2_b = 3*pi*ρ*b^2

    d = b/2 - a*b
    d_b = 1/2 - a

    L = tmp1*(v + d*ω - u*α0) + tmp2*u/b*ω

    L_b = tmp1_b*(v + d*ω - u*α0) + tmp1*d_b*ω + tmp2_b*u/b*ω - tmp2*u/b^2*ω
    M_b = -tmp2_b*u*ω + (1/2 + a)*L + (b/2 + a*b)*L_b

    return SVector(L_b, M_b)
end

function quasisteady1_loads_ρ(a, b, ρ, a0, α0, u, v, ω)

    tmp1_ρ = a0*u*b

    tmp2_ρ = pi*b^3

    d = b/2 - a*b

    L_ρ = tmp1_ρ*(v + d*ω - u*α0) + tmp2_ρ*u/b*ω

    M_ρ = -tmp2_ρ*u*ω + (b/2 + a*b)*L_ρ

    return SVector(L_ρ, M_ρ)
end

function quasisteady1_loads_a0(a, b, ρ, a0, α0, u, v, ω)

    tmp1_a0 = ρ*u*b

    tmp2 = pi*ρ*b^3

    d = b/2 - a*b

    L_a0 = tmp1_a0*(v + d*ω - u*α0)
    M_a0 = (b/2 + a*b)*L_a0

    return SVector(L_a0, M_a0)
end

function quasisteady1_loads_α0(a, b, ρ, a0, u)

    tmp1 = a0*ρ*u*b

    L_α0 = -tmp1*u

    M_α0 = (b/2 + a*b)*L_α0

    return SVector(L_α0, M_α0)
end

function quasisteady1_u(a, b, ρ, a0, α0, u, v, ω)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    tmp1_u = a0*ρ*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L_u = tmp1_u*(v + d*ω - u*α0) - tmp1*α0 + tmp2/b*ω
    # moment at reference point
    M_u = -tmp2*ω + (b/2 + a*b)*L_u

    return SVector(L_u, M_u)
end

function quasisteady1_v(a, b, ρ, a0, α0, u)
    # lift at reference point
    L_v = a0*ρ*u*b
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

function quasisteady1_ω(a, b, ρ, a0, u)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L_ω = tmp1*d + tmp2*u/b
    # moment at reference point
    M_ω = -tmp2*u + (b/2 + a*b)*L_ω

    return SVector(L_ω, M_ω)
end

# quasi-steady loads + added mass effects
function quasisteady2_loads(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L = tmp1*(v + d*ω - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + (b/8 - a*b/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function quasisteady2_loads_a(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)

    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3

    L = tmp1*(v + (b/2 - a*b)*ω - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    L_a = -tmp1*b*ω - tmp2*ωdot

    M_a = tmp2*b/2*ωdot + b*L + (b/2 + a*b)*L_a

    return SVector(L_a, M_a)
end

function quasisteady2_loads_b(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)

    tmp1 = a0*ρ*u*b
    tmp1_b = a0*ρ*u

    tmp2 = pi*ρ*b^3
    tmp2_b = 3*pi*ρ*b^2

    d = b/2 - a*b
    d_b = 1/2 - a

    L = tmp1*(v + d*ω - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)

    L_b = tmp1_b*(v + d*ω - u*α0) + tmp1*d_b*ω +
        tmp2_b*(vdot/b + u/b*ω - a*ωdot) + tmp2*(-vdot/b^2 - u/b^2*ω)
    M_b = -tmp2_b*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) - tmp2*(1/8 - a/2)*ωdot +
        (1/2 + a)*L + (b/2 + a*b)*L_b

    return SVector(L_b, M_b)
end

function quasisteady2_loads_ρ(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)

    tmp1_ρ = a0*u*b

    tmp2_ρ = pi*b^3

    d = b/2 - a*b

    L_ρ = tmp1_ρ*(v + d*ω - u*α0) + tmp2_ρ*(vdot/b + u/b*ω - a*ωdot)

    M_ρ = -tmp2_ρ*(vdot/2 + u*ω + (b/8 - a*b/2)*ωdot) + (b/2 + a*b)*L_ρ

    return SVector(L_ρ, M_ρ)
end

quasisteady2_loads_a0(args...) = quasisteady1_loads_a0(args...)

quasisteady2_loads_α0(args...) = quasisteady1_loads_α0(args...)

quasisteady2_loads_u(args...) = quasisteady1_u(args...)

quasisteady2_loads_v(args...) = quasisteady1_v(args...)

quasisteady2_loads_ω(args...) = quasisteady1_ω(args...)

function quasisteady2_loads_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_vdot = tmp/b
    # moment at reference point
    M_vdot = -tmp/2 + (b/2 + a*b)*L_vdot

    return SVector(L_vdot, M_vdot)
end

function quasisteady2_loads_ωdot(a, b, ρ)
    tmp = pi*ρ*b^3
    L_ωdot = -tmp*a
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_ωdot
    return SVector(L_ωdot, M_ωdot)
end

function quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, ω)

    return quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
end

# quasi-steady loads from acceleration terms
function quasisteady2_rate_loads(a, b, ρ, vdot, ωdot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(vdot/b - a*ωdot)
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end
