"""
    PetersFiniteState{N, TF} <: AerodynamicModel

Peter's finite state model with `N` aerodynamic states.
"""
struct PetersFiniteState{N, TF} <: StructuralModel
    A::SMatrix{N,N,TF}
    b::SVector{N,TF}
    c::SVector{N,TF}
end

"""
    PetersFiniteState{N,TF=Float64}()

Initialize an object of type `PetersFiniteState` which has `N` aerodynamic
degrees of freedom.
"""
PetersFiniteState{N}() where N = PetersFiniteState{N,Float64}()

function PetersFiniteState{N,TF}() where {N, TF}

    b = zeros(TF, N)
    for i = 1:N
        b[i] = (-1)^(i-1)*factorial(big(N + i))/factorial(big(N - i))*
            1/factorial(big(i))^2
    end
    b[N] += (-1)^N

    c = zeros(TF, N)
    c .= 2/N

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[m, n] = n/2
    end
    for m in 2:N
        n = m - 1
        D[m, n] = -n/2
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return PetersFiniteState{N, TF}(A, b, c)
end

isinplace(::PetersFiniteState, ::TypicalSection) = false
identity_mass_matrix(::PetersFiniteState, ::TypicalSection) = false
constant_aero_mass_matrices(::PetersFiniteState, ::TypicalSection) = false
constant_load_mass_matrices(::PetersFiniteState, ::TypicalSection) = false
defined_state_jacobians(::PetersFiniteState, ::TypicalSection) = true
defined_load_jacobians(::PetersFiniteState, ::TypicalSection) = true

function init_aero_mass_matrices(::PetersFiniteState{N,TF}, ::TypicalSection) where {N, TF}
    Mas = zeros(TF, N, 4)
    Maa = zeros(TF, N, N)
    return Mas, Maa
end

function update_aero_mass_matrices!(aero::PetersFiniteState, stru::TypicalSection,
    Mas, Maa, q, λ, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    # extract model constants
    Abar, cbar = aero.A, aero.c
    # update mass matrices
    Mas .= peters_stru_mass_matrix(a, b, cbar)
    Maa .= peters_aero_mass_matrix(Abar)
    # return updated mass matrices
    return Mas, Maa
end

function get_aero_rates(aero::PetersFiniteState, stru::TypicalSection, q, λ, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    # extract model constants
    cbar = aero.c
    # calculate rates
    return peters_rhs(b, U, cbar, θdot, λ)
end

function get_aero_jacobians(aero::PetersFiniteState{N,TF}, stru::TypicalSection, q, λ,
    p, t) where {N, TF}
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    # extract model parameters
    cbar = aero.c
    # jacobian with respect to structural states
    Jas = peters_stru_jacobian(U, cbar)
    # jacobian with respect to aerodynamic states
    Jaa = peters_aero_jacobian(b, U, cbar)
    # return jacobians
    return Jas, Jaa
end

function init_load_mass_matrices(::PetersFiniteState{N,TF}, ::TypicalSection) where {N, TF}
    Mls = zeros(TF, 2, 4)
    Mla = zeros(TF, 2, N)
    return Mas, Maa
end

function update_load_mass_matrices!(aero::PetersFiniteState{N,TF},
    stru::TypicalSection, Mls, Mla, q, λ, p, t) where {N,TF}
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    # extract model parameters
    cbar = aero.c
    # update mass matrices
    Mls .= peters_load_stru_mass_matrix(a, b, ρ)
    Mla .= peters_load_aero_mass_matrix(cbar)
    # return updated mass matrices
    return Mls, Mla
end

function get_loads(aero::PetersFiniteState, stru::TypicalSection, q, λ, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    # extract model parameters
    bbar = aero.b
    # calculate aerodynamic loads
    return peters_state_loads(a, b, U, ρ, bbar, θ, hdot, θdot, λ)
end

function get_load_jacobians(aero::PetersFiniteState, stru::TypicalSection, q, λ, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # extract model parameters
    bbar = aero.b
    # jacobian with respect to structural states
    Jls  = peters_load_stru_jacobian(a, b, U, ρ)
    # jacobian with respect to aerodynamic states
    Jla = peters_load_aero_jacobian(b, U, ρ, bbar)
    # return jacobians
    return Jls, Jla
end

# --- Internal Functions --- #
peters_stru_mass_matrix(a, b, cbar) = hcat(zero(cbar), zero(cbar), -cbar, -cbar*(b/2-a*b))
peters_aero_mass_matrix(Abar) = Abar
peters_lhs(a, b, Abar, cbar, dhdot, dθdot, dλ) = Abar*dλ - cbar*(dhdot + b*(1/2 - a)*dθdot)
peters_rhs(b, U, cbar, θdot, λ) = cbar*U*θdot - U/b*λ
peters_stru_jacobian(U, cbar) = hcat(zero(cbar), zero(cbar), zero(cbar), U*cbar)
peters_aero_jacobian(b, U, cbar) = -U/b*Diagonal(one.(cbar))
peters_load_stru_mass_matrix(a, b, ρ) = pi*ρ*b^2*@SMatrix [0 0 -1 b*a; 0 0 b/2 b^2*(1/8 - a/2)]
peters_load_aero_mass_matrix(cbar) = vcat(zero(cbar)', zero(cbar)')
function peters_state_loads(a, b, U, ρ, bbar, θ, hdot, θdot, λ)
    # calculate induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # calculate lift (excluding state rate terms)
    L = 2*pi*ρ*U*b*(hdot + U*θ + (b/2-a*b)*θdot) + pi*ρ*b^2*(U*θdot - λ0)
    # calculate moment (excluding state rate terms)
    M = -pi*ρ*b^3*U*θdot
    # return load
    return SVector(L, M)
end
function peters_rate_loads(a, b, ρ, dhdot, dθdot)
    L = pi*ρ*b^2*(dhdot - b*a*dθdot)
    M = -pi*ρ*b^3*(1/2*dhdot + b*(1/8-a/2)*dθdot)
    return SVector(L, M)
end
peters_load_stru_jacobian(a, b, U, ρ) = 2*pi*ρ*b*U*(@SMatrix [0 U 1 b-a*b; 0 0 0 -b^2/2])
peters_load_aero_jacobian(b, ρ, bbar) = -pi/2*ρ*b^2*vcat(bbar', zero(bbar)')
