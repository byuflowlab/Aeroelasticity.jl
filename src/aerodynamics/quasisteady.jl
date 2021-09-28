"""
    QuasiSteady{Order} <: NoStateModel

2D quasi-steady aerodynamic model with parameters ``a, b, a_0, \\alpha_0, c_{d0}, c_{m0}``.
"""
struct QuasiSteady{Order} <: NoStateModel end

"""
    Steady <: NoStateModel

2D steady aerodynamic model with parameters ``a, b, a_0, \\alpha_0, c_{d0}, c_{m0}``.
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

number_of_parameters(::Type{<:QuasiSteady}) = 6

# --- Convenience Methods --- #

function set_parameters!(p, model::QuasiSteady; a, b, a0, alpha0, cd0, cm0)

    p[1] = a
    p[2] = b
    p[3] = a0
    p[4] = alpha0
    p[5] = cd0
    p[6] = cm0

    return p
end

function separate_parameters(model::QuasiSteady, p)

    return (a=p[1], b=p[2], a0=p[3], alpha0=p[4], cd0=p[5], cm0=p[6])
end

# --- Internal Methods --- #

# steady state loads
function quasisteady0_loads(a, b, ρ, a0, α0, cd0, cm0, u, v)
    # normal force at reference point
    N = a0*ρ*b*u*(v - u*α0)
    # axial force at reference point
    A = -a0*ρ*b*(v - u*α0)^2 + ρ*b*u^2*cd0
    # moment at reference point
    M = 2*ρ*b^2*u^2*cm0 + (b/2 + a*b)*N

    return SVector(N, A, M)
end

function quasisteady0_loads_a(a, b, ρ, a0, α0, u, v)

    N = a0*ρ*b*u*(v - u*α0)

    N_a = 0

    A_a = 0

    M_a = b*N

    return SVector(N_a, A_a, M_a)
end

function quasisteady0_loads_b(a, b, ρ, a0, α0, cd0, cm0, u, v)

    N = a0*ρ*b*u*(v - u*α0)

    N_b = a0*ρ*u*(v - u*α0)

    A_b = -a0*ρ*(v - u*α0)^2 + ρ*u^2*cd0

    M_b = 4*ρ*b*u^2*cm0 + (1/2 + a)*N + (b/2 + a*b)*N_b

    return SVector(N_b, A_b, M_b)
end

function quasisteady0_loads_ρ(a, b, ρ, a0, α0, cd0, cm0, u, v)
    N_ρ = a0*b*u*(v - u*α0)
    A_ρ = -a0*b*(v - u*α0)^2 + b*u^2*cd0
    M_ρ = 2*b^2*u^2*cm0 + (b/2 + a*b)*N_ρ
    return SVector(N_ρ, A_ρ, M_ρ)
end

function quasisteady0_loads_a0(a, b, ρ, α0, u, v)
    # normal force at reference point
    N_a0 = ρ*b*u*(v - u*α0)
    A_a0 = -ρ*b*(v - u*α0)^2
    M_a0 = (b/2 + a*b)*N_a0

    return SVector(N_a0, A_a0, M_a0)
end

function quasisteady0_loads_α0(a, b, ρ, a0, α0, u, v)
    N_α0 = -a0*ρ*b*u^2
    A_α0 = 2*a0*b*u*(v - u*α0)
    M_α0 = (b/2 + a*b)*N_α0

    return SVector(N_α0, A_α0, M_α0)
end

function quasisteady0_loads_cd0(b, ρ, u)
    N_cd0 = 0
    A_cd0 = ρ*b*u^2
    M_cd0 = 0
    return SVector(N_cd0, A_cd0, M_cd0)
end

function quasisteady0_loads_cm0(b, ρ, u)
    N_cm0 = 0
    A_cm0 = 0
    M_cm0 = 2*ρ*b^2*u^2
    return SVector(N_cm0, A_cm0, M_cm0)
end

function quasisteady0_loads_u(a, b, ρ, a0, α0, cd0, cm0, u, v)
    N_u = a0*ρ*b*v - 2*a0*ρ*b*u*α0
    A_u = 2*a0*α0*ρ*b*(v - u*α0)^2 + 2*ρ*b*u*cd0
    M_u = 4*ρ*b^2*u*cm0 + (b/2 + a*b)*N_u
    return SVector(N_u, A_u, M_u)
end

function quasisteady0_loads_v(a, b, ρ, a0, α0, u, v)
    N_v = a0*ρ*b*u
    A_v = -2*a0*ρ*b*(v - u*α0)
    M_v = (b/2 + a*b)*N_v
    return SVector(N_v, A_v, M_v)
end

# circulatory loads + added mass effects, neglecting accelerations
function quasisteady1_loads(a, b, ρ, a0, α0, cd0, cm0, u, v, ω)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # normal force at reference point
    N = tmp1*(v + d*ω - u*α0) + tmp2*u/b*ω
    # axial force at the reference point
    A = -a0*ρ*b*(v + d*ω - u*α0)^2 + ρ*b*u^2*cd0
    # moment at reference point
    M = -tmp2*u*ω + 2*ρ*b^2*u^2*cm0 + (b/2 + a*b)*N

    return SVector(N, A, M)
end

function quasisteady1_loads_a(a, b, ρ, a0, α0, u, v, ω)

    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3

    N = tmp1*(v + (b/2 - a*b)*ω - u*α0) + tmp2*u/b*ω
    N_a = -tmp1*b*ω
    A_a = 2*a0*ρ*b^2*(v + (b/2 - a*b)*ω - u*α0)*ω
    M_a = b*N + (b/2 + a*b)*N_a

    return SVector(N_a, A_a, M_a)
end


function quasisteady1_loads_b(a, b, ρ, a0, α0, cd0, cm0, u, v, ω)

    tmp1 = a0*ρ*u*b
    tmp1_b = a0*ρ*u

    tmp2 = pi*ρ*b^3
    tmp2_b = 3*pi*ρ*b^2

    d = b/2 - a*b
    d_b = 1/2 - a

    N = tmp1*(v + d*ω - u*α0) + tmp2*u/b*ω

    N_b = tmp1_b*(v + d*ω - u*α0) + tmp1*d_b*ω + tmp2_b*u/b*ω - tmp2*u/b^2*ω
    A_b = -2*a0*ρ*b*(v + d*ω - u*α0)*d_b - a0*ρ*(v + d*ω - u*α0)^2  + ρ*u^2*cd0
    M_b = -tmp2_b*u*ω + 4*ρ*b*u^2*cm0 + (1/2 + a)*N + (b/2 + a*b)*N_b

    return SVector(N_b, A_b, M_b)
end

function quasisteady1_loads_ρ(a, b, ρ, a0, α0, cd0, cm0, u, v, ω)

    tmp1_ρ = a0*u*b

    tmp2_ρ = pi*b^3

    d = b/2 - a*b

    N_ρ = tmp1_ρ*(v + d*ω - u*α0) + tmp2_ρ*u/b*ω
    A_ρ = -a0*b*(v + d*ω - u*α0)^2 + b*u^2*cd0
    M_ρ = -tmp2_ρ*u*ω + 2*b^2*u^2*cm0 + (b/2 + a*b)*N_ρ

    return SVector(N_ρ, A_ρ, M_ρ)
end

function quasisteady1_loads_a0(a, b, ρ, a0, α0, u, v, ω)

    tmp1_a0 = ρ*u*b

    tmp2 = pi*ρ*b^3

    d = b/2 - a*b

    N_a0 = tmp1_a0*(v + d*ω - u*α0)
    A_a0 = -ρ*b*(v + d*ω - u*α0)^2
    M_a0 = (b/2 + a*b)*N_a0

    return SVector(N_a0, A_a0, M_a0)
end

function quasisteady1_loads_α0(a, b, ρ, a0, α0, u, v, ω)
    tmp1 = a0*ρ*u*b
    N_α0 = -tmp1*u
    A_α0 = 2*a0*ρ*b*u*(v + (b/2 - a*b)*ω - u*α0)
    M_α0 = (b/2 + a*b)*N_α0
    return SVector(N_α0, A_α0, M_α0)
end

function quasisteady1_loads_cd0(b, ρ, u)
    N_cd0 = 0
    A_cd0 = ρ*b*u^2
    M_cd0 = 0
    return SVector(N_cd0, A_cd0, M_cd0)
end

function quasisteady1_loads_cm0(b, ρ, u)
    N_cm0 = 0
    A_cm0 = 0
    M_cm0 = 2*ρ*b^2*u^2
    return SVector(N_cm0, A_cm0, M_cm0)
end

function quasisteady1_loads_u(a, b, ρ, a0, α0, cd0, cm0, u, v, ω)
    tmp1 = a0*ρ*u*b
    tmp1_u = a0*ρ*b
    tmp2 = pi*ρ*b^3
    d = b/2 - a*b
    N_u = tmp1_u*(v + d*ω - u*α0) - tmp1*α0 + tmp2/b*ω
    A_u = 2*a0*ρ*b*α0*(v + d*ω - u*α0) + 2*ρ*b*u*cd0
    M_u = -tmp2*ω + (b/2 + a*b)*N_u + 4*ρ*b^2*u*cm0
    return SVector(N_u, A_u, M_u)
end

function quasisteady1_loads_v(a, b, ρ, a0, α0, u, v, ω)
    d = b/2 - a*b
    N_v = a0*ρ*u*b
    A_v = -2*a0*ρ*b*(v + d*ω - u*α0)
    M_v = (b/2 + a*b)*N_v
    return SVector(N_v, A_v, M_v)
end

function quasisteady1_loads_ω(a, b, ρ, a0, u)
    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3
    d = b/2 - a*b
    N_ω = tmp1*d + tmp2*u/b
    A_ω = 2*a0*ρ*a*b^2*(v + d*ω - u*α0)
    M_ω = -tmp2*u + (b/2 + a*b)*N_ω
    return SVector(N_ω, A_ω, M_ω)
end

# quasi-steady loads + added mass effects
function quasisteady2_loads(a, b, ρ, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # normal force at reference point
    N = tmp1*(v + d*ω - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # axial force at reference point
    A = -a0*ρ*b*(v + d*ω - u*α0)^2 + ρ*b*u^2*cd0
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + (b/8 - a*b/2)*ωdot) + 2*ρ*b^2*u^2*cm0 + (b/2 + a*b)*N

    return SVector(N, A, M)
end

function quasisteady2_loads_a(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)

    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3

    N = tmp1*(v + (b/2 - a*b)*ω - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    N_a = -tmp1*b*ω - tmp2*ωdot
    A_a = 2*a0*ρ*b^2*(v + (b/2 - a*b)*ω - u*α0)*ω
    M_a = tmp2*b/2*ωdot + b*N + (b/2 + a*b)*N_a

    return SVector(N_a, A_a, M_a)
end

function quasisteady2_loads_b(a, b, ρ, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)

    tmp1 = a0*ρ*u*b
    tmp1_b = a0*ρ*u

    tmp2 = pi*ρ*b^3
    tmp2_b = 3*pi*ρ*b^2

    d = b/2 - a*b
    d_b = 1/2 - a

    N = tmp1*(v + d*ω - u*α0) + tmp2*(vdot/b + u/b*ω - a*ωdot)

    N_b = tmp1_b*(v + d*ω - u*α0) + tmp1*d_b*ω +
        tmp2_b*(vdot/b + u/b*ω - a*ωdot) + tmp2*(-vdot/b^2 - u/b^2*ω)
    A_b = -2*a0*ρ*b*(v + d*ω - u*α0)*d_b - a0*ρ*(v + d*ω - u*α0)^2 + ρ*u^2*cd0
    M_b = -tmp2_b*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) - tmp2*(1/8 - a/2)*ωdot +
        4*ρ*b*u^2*cm0 + (1/2 + a)*N + (b/2 + a*b)*N_b

    return SVector(N_b, A_b, M_b)
end

function quasisteady2_loads_ρ(a, b, ρ, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)

    tmp1_ρ = a0*u*b

    tmp2_ρ = pi*b^3

    d = b/2 - a*b

    N_ρ = tmp1_ρ*(v + d*ω - u*α0) + tmp2_ρ*(vdot/b + u/b*ω - a*ωdot)
    A_ρ = -a0*b*(v + d*ω - u*α0)^2 + b*u^2*cd0
    M_ρ = -tmp2_ρ*(vdot/2 + u*ω + (b/8 - a*b/2)*ωdot) + 2*b^2*u^2*cm0 + (b/2 + a*b)*N_ρ

    return SVector(N_ρ, A_ρ, M_ρ)
end

quasisteady2_loads_a0(args...) = quasisteady1_loads_a0(args...)

quasisteady2_loads_α0(args...) = quasisteady1_loads_α0(args...)

quasisteady2_loads_cd0(args...) = quasisteady1_loads_cd0(args...)

quasisteady2_loads_cm0(args...) = quasisteady1_loads_cm0(args...)

quasisteady2_loads_u(args...) = quasisteady1_loads_u(args...)

quasisteady2_loads_v(args...) = quasisteady1_loads_v(args...)

quasisteady2_loads_ω(args...) = quasisteady1_loads_ω(args...)

function quasisteady2_loads_vdot(a, b, ρ)
    tmp = pi*ρ*b^3
    N_vdot = tmp/b
    A_vdot = 0
    M_vdot = -tmp/2 + (b/2 + a*b)*N_vdot
    return SVector(N_vdot, A_vdot, M_vdot)
end

function quasisteady2_loads_ωdot(a, b, ρ)
    tmp = pi*ρ*b^3
    N_ωdot = -tmp*a
    A_ωdot = 0
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*N_ωdot
    return SVector(N_ωdot, A_ωdot, M_ωdot)
end