"""
    Wagner{N, TF} <: AerodynamicModel

Aerodynamic model based on Wagner's function, with two degrees of freedom.
"""
struct Wagner <: AerodynamicModel end
isinplace(::Wagner) = false
islinear(::Wagner) = true
is2D(::Wagner) = true

# --- Equations of Motion and Jacobians --- #

wagner_lhs(dλ) = dλ

function wagner_rhs(a, b, U, C1, C2, ε1, ε2, hdot, θdot, λ1, λ2)
    λ1dot = -ε1*U/b*λ1 + C1*ε1*U/b*(hdot + (1/2-a)*b*θdot + U*θ)
    λ2dot = -ε2*U/b*λ2 + C2*ε2*U/b*(hdot + (1/2-a)*b*θdot + U*θ)
    return SVector(λ1dot, λ2dot)
end

wagner_rates(args...) = wagner_rhs(args...)

wagner_jacobian(a, b, U, ε1, ε2) = @SMatrix [-ε1*V/b 0; 0 -ε2*V/b]

wagner_mass_matrix() = I

wagner_coupled_lhs(dλ) = dλ

function wagner_coupled_rhs(a, b, U, C1, C2, ε1, ε2, hdot, θdot, λ1, λ2)
    return wagner_rhs(a, b, U, C1, C2, ε1, ε2, hdot, θdot, λ1, λ2)
end

wagner_aerodynamic_jacobian(a, b, U, ε1, ε2) = wagner_jacobian(a, b, U, ε1, ε2)

wagner_aerodynamic_mass_matrix() = I

function wagner_structural_jacobian(a, b, U, C1, C2, ε1, ε2)
    tmp1 = C1*ε1*U/b
    tmp2 = C2*ε2*U/b
    tmp3 = (1/2-a)*b
    return @SMatrix [
        0 0 0 0;
        0 0 0 0;
        0 tmp1*V tmp1 tmp1*tmp3;
        0 tmp2*V tmp2 tmp2*tmp3;
        ]
end

wagner_structural_mass_matrix() = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]

function wagner_loads(a, b, U, ρ, C1, C2, θ, hdot, θdot, hddot, θddot, λ)
    ϕ0 = 1 - C1 - C2
    wt = hdot + (1/2 - a)*b*θdot + U*θ
    L = pi*ρ*b^2*(hddot + U*θdot - a*b*θddot) +
        2*pi*ρ*b*U*(wt*ϕ0 + λ1 + λ2)
    M = pi*ρ*b^2*(a*b*hddot - (1/2-a)*b*U*θdot - (1/8 + a^2)*b^2*θddot) +
        2*pi*ρ*b^2*U*(1/2+a)*(wt*ϕ0 + λ1 + λ2)
    return SVector(L, M)
end

function peters_loads_lhs(a, b, ρ, dhdot, dθdot)
    SVector(
        pi*ρ*b^2*(b*a*dθdot - dhdot), # lift moved to lhs
        pi*ρ*b^3*(1/2*dhdot + b*(1/8 - a/2)*dθdot) # moment moved to lhs
    )
end

function peters_loads_rhs(a, b, U, ρ, bn, θ, hdot, θdot, λ)
    λ0 = 1/2 * bn'*λ # induced flow velocity
    SVector(
        2*pi*ρ*U*b*(hdot + U*θ + b*(1/2-a)*θdot) + pi*ρ*b^2*(U*θdot - λ0), # lift
        -pi*ρ*b^3*U*θdot # moment
    )
end

# jacobian of loads wrt aerodynamic states
function peters_loads_aerodynamic_jacobian(b, ρ, bn)
    return vcat(-1/2*pi*ρ*b^2*bn', zero(bn)')
end

peters_loads_aerodynamic_mass_matrix(bn) = vcat(zero(bn)', zero(bn)')

function peters_loads_structural_jacobian(a, b, U, ρ)
    tmp1 = 2*pi*ρ*U*b
    tmp2 = pi*ρ*b^2
    tmp3 = -pi*ρ*b^3
    return @SMatrix [0 tmp1*U tmp1 tmp1*b*(1/2-a)+tmp2*U; 0 0 0 tmp3*U]
end

function peters_loads_structural_mass_matrix(a, b, U, ρ)
    tmp1 = 2*pi*ρ*U*b
    tmp2 = pi*ρ*b^2
    tmp3 = -pi*ρ*b^3
    return @SMatrix [0 0 -tmp2 tmp2*b*a; 0 0 -tmp3*1/2 -tmp3*b*(1/8-a/2)]
end

# mass matrix multiplied loads
function peters_coupled_loads(U, ρ, a, b, bn, θ, hdot, θdot, λ)
    λ0 = 1/2 * bn'*λ # induced flow velocity
    L = 2*pi*ρ*U*b*(hdot + U*θ + b*(1/2-a)*θdot) + # circulatory lift
        pi*ρ*b^2*(U*θdot - λ0) # non-circulatory lift
    M = -pi*ρ*b^3*U*θdot # moment
    return SVector(L, M)
end

# --- Aerodynamic Interface --- #

"""
    get_ode(::Wagner)

Return an ODE function of the form `du = f(u, p, t)` for Peter's finite state
model with states `u = [λ1, λ2, ... λN]` and parameters `p = [U, a, b, c, θdot(t),
hddot(t), θddot(t)]`
"""
get_ode(model::Wagner) = (u, p, t) -> get_aerodynamic_rates(model, u, p, t)

function get_aerodynamic_rates(model::Wagner, u, p, t)
    λ = u
    a, b, U, θdot, hddot, θddot = p
    Ainv, c = model.Ainv, model.c
    θdot = typeof(θdot) <: Number ? θdot : θdot(t)
    hddot = typeof(hddot) <: Number ? hddot : hddot(t)
    θddot = typeof(θddot) <: Number ? θddot : θddot(t)
    return peters_rates(a, b, U, Ainv, c, θdot, hddot, θddot, λ)
end

get_mass_matrix(::Wagner) = I

function get_jacobian(model::Wagner, u, p, t)
    a, b, U = p
    Ainv = model.Ainv
    return peters_jacobian(b, U, Ainv)
end

# TODO: Add get_parameter_jacobian

@generated function get_syms(::Wagner{N, <:Any}) where N
    syms = Tuple(Symbol(:λ, i) for i = 1:N)
    return :($syms)
end

# --- TypicalSection Interface --- #

# mass matrix multiplied rates
function get_aerodynamic_rhs(aero::Wagner{N, <:Any}, stru::TypicalSection,
    u, p, t) where N

    # extract states
    iu_a = SVector{4}(1:4)
    iu_λ = SVector{N}(5:4+N)
    h, θ, hdot, θdot = u[iu_a]
    λ = u[iu_λ]

    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p

    # return mass matrix multiplied rates
    return peters_coupled_rhs(a, b, U, aero.c, θdot, λ)
end

function get_loads(aero::Wagner{N, <:Any}, stru, u, p, t) where N

    # extract states
    iu_a = SVector{4}(1:4)
    iu_λ = SVector{N}(5:4+N)
    h, θ, hdot, θdot = u[iu_a]
    λ = u[iu_λ]

    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p

    # return loads
    return peters_coupled_loads(U, ρ, a, b, aero.b, θ, hdot, θdot, λ)
end

function get_load_mass_matrices(aero::Wagner, stru::TypicalSection, u, p, t)
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    Mls = peters_loads_structural_mass_matrix(a, b, U, ρ)
    Mla = peters_loads_aerodynamic_mass_matrix(aero.b)
    return Mls, Mla
end

function get_aerodynamic_mass_matrices(aero::Wagner, stru::TypicalSection, u, p, t)
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    Mas = peters_structural_mass_matrix(a, b, aero.c)
    Maa = peters_aerodynamic_mass_matrix(aero.A)
    return Mas, Maa
end

function get_load_jacobians(aero::Wagner, stru::TypicalSection, u, p, t)
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    Jls = peters_loads_structural_jacobian(a, b, U, ρ)
    Jla = peters_loads_aerodynamic_jacobian(b, ρ, aero.b)
    return Jls, Jla
end

function get_aerodynamic_jacobians(aero::Wagner, stru::TypicalSection, u, p, t)
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    Jas = peters_structural_jacobian(U, aero.c)
    Jaa = peters_aerodynamic_jacobian(b, U, aero.c)
    return Jas, Jaa
end
