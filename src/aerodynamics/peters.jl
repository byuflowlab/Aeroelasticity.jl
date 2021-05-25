"""
    PetersFiniteState{N,TF} <: AerodynamicModel

Peter's finite state model with `N` state variables, inputs ``d = \\begin{bmatrix}
\\dot{\\theta} & \\ddot{h} & \\ddot{\\theta}\\end{bmatrix}^T`` and parameters
``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``
"""
struct PetersFiniteState{N,TF} <: StructuralModel
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

function PetersFiniteState{N,TF}() where {N,TF}

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

    return PetersFiniteState{N,TF}(A, b, c)
end

number_of_states(::PetersFiniteState{N,TF}) where {N,TF} = N
number_of_inputs(::PetersFiniteState) = 3
number_of_parameters(::PetersFiniteState) = 4
isinplace(::PetersFiniteState) = false
has_mass_matrix(::PetersFiniteState) = true
constant_mass_matrix(::PetersFiniteState) = true
linear_input_dependence(::PetersFiniteState) = true
defined_state_jacobian(::PetersFiniteState) = true
defined_input_jacobian(::PetersFiniteState) = true
constant_input_jacobian(::PetersFiniteState) = false

get_mass_matrix(model::PetersFiniteState) =  model.A

function get_rates(model::PetersFiniteState, λ, d, p, t)
    # extract structural deflections
    θdot, hddot, θddot = d
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return peters_rhs(b, U, cbar, θdot, λ)
end

function get_state_jacobian(model::PetersFiniteState, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return peters_state_jacobian(b, U, cbar)
end

function get_input_jacobian(model::PetersFiniteState, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return peters_input_jacobian(a, b, U, cbar)
end

# TODO: Add parameter jacobian

# --- Internal Functions --- #
peters_lhs(a, b, Abar, cbar, dhdot, dθdot, dλ) = Abar*dλ
peters_rhs(b, U, cbar, θdot, λ) = cbar*(hddot + U*θdot + (b/2-a*b)*θddot) - U/b*λ
peters_mass_matrix(Abar) = Abar
peters_state_jacobian(b, U, cbar) = -U/b*Diagonal(one.(cbar))
peters_input_jacobian(a, b, U, cbar) = hcat(U*cbar, cbar, (b/2-a*b)*cbar)

function peters_loads(a, b, U, ρ, bbar, θ, hdot, θdot, λ)
    # calculate induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # calculate lift (excluding state rate terms)
    L = 2*pi*ρ*U*b*(hdot + U*θ + (b/2-a*b)*θdot) + pi*ρ*b^2*(U*θdot - λ0)
    # calculate moment (excluding state rate terms)
    M = -pi*ρ*b^3*U*θdot
    # return load
    return SVector(L, M)
end
peters_loads_λ(b, ρ, bbar) = -pi/2*ρ*b^2*vcat(bbar', zero(bbar)')
function peters_loads_λdot(cbar)
    tmp = zero(cbar)'
    return vcat(tmp, tmp)
end
peters_loads_h(b, U, ρ) = SVector(2*pi*ρ*U*b, 0)
peters_loads_θ(b, U, ρ) = SVector(2*pi*ρ*U^2*b, 0)
peters_loads_hdot() = SVector(0, 0)
peters_loads_θdot(a, b, U, ρ) = SVector(pi/2*ρ*b^2*U*(3 - 4*a), -pi*ρ*b^3*U)
peters_loads_hddot(b, ρ) = SVector(pi*ρ*b^2, -pi/2*ρ*b^3)
peters_loads_θddot(a, b, ρ) = SVector(-pi*ρ*a*b^3, -pi/8*ρ*b^4*(1 - 4*a))
