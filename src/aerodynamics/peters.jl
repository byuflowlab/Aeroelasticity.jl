"""
    PetersFiniteState{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``d = \\begin{bmatrix}
\\dot{\\theta} & \\ddot{h} & \\ddot{\\theta}\\end{bmatrix}^T`` and parameters
``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``
"""
struct PetersFiniteState{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

# --- Constructors --- #
"""
    PetersFiniteState{N,TF=Float64}()

Initialize an object of type `PetersFiniteState` which has `N` aerodynamic
degrees of freedom.
"""
PetersFiniteState{N}() where N = PetersFiniteState{N,Float64}()

function PetersFiniteState{N,TF}() where {N,TF}

    b = zeros(TF, N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(N + n - 1)/factorial(N - n - 1)*1/factorial(n)^2
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
        D[m, n] = 1/(2*n)
    end
    for m in 2:N
        n = m - 1
        D[m, n] = -1/(2*n)
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return PetersFiniteState(SMatrix{N,N,TF}(A), SVector{N,TF}(b), SVector{N,TF}(c))
end

# --- Traits --- #
number_of_states(::Type{PetersFiniteState{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
number_of_inputs(::Type{<:PetersFiniteState}) = 3
number_of_parameters(::Type{<:PetersFiniteState}) = 4
isinplace(::Type{<:PetersFiniteState}) = false
has_mass_matrix(::Type{<:PetersFiniteState}) = true
constant_mass_matrix(::Type{<:PetersFiniteState}) = true
linear_input_dependence(::Type{<:PetersFiniteState}) = true
defined_state_jacobian(::Type{<:PetersFiniteState}) = true
defined_input_jacobian(::Type{<:PetersFiniteState}) = true
constant_input_jacobian(::Type{<:PetersFiniteState}) = false

# --- Interface Methods --- #
get_mass_matrix(model::PetersFiniteState) = model.A

function get_rates(model::PetersFiniteState{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract structural deflections
    θdot, hddot, θddot = d
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return peters_rhs(a, b, U, cbar, θdot, hddot, θddot, λ)
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

# --- Typical Section Model Coupling --- #

# traits
inplace_inputs(::Type{<:PetersFiniteState}, ::Type{TypicalSection}) = false
has_input_mass_matrix(::Type{<:PetersFiniteState}, ::Type{TypicalSection}) = true
constant_input_mass_matrix(::Type{<:PetersFiniteState}, ::Type{TypicalSection}) = false
defined_input_state_jacobian(::Type{<:PetersFiniteState}, ::Type{TypicalSection}) = true

# interface methods
function get_input_mass_matrix(aero::PetersFiniteState, stru::TypicalSection, u, p, t)
    # extract model constants
    cbar = aero.c
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # create zero row vector with length nλ
    zλ = zero(aero.c')
    # construct submatrices
    Mds = @SMatrix [0 0 0 0; 0 0 -1 0; 0 0 0 -1]
    Mda = vcat(zλ, zλ, zλ)
    Mrs = hcat(
        (@SVector zeros(2)),
        (@SVector zeros(2)),
        peters_loads_hddot(b, ρ),
        peters_loads_θddot(a, b, ρ),
        )
    Mra = peters_loads_λdot(aero.c)
    # assemble mass matrix
    return vcat(hcat(Mds, Mda), hcat(Mrs, Mra))
end

function get_inputs(aero::PetersFiniteState{N,TF,SV,SA}, stru::TypicalSection,
    u, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iq = SVector{4}(1:4)
    iλ = SVector{N}(5:4+N)
    # separate aerodynamic and structural states
    q = u[iq]
    λ = u[iλ]
    # extract structural state variables
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # extract model constants
    bbar = aero.b
    # calculate aerodynamic loads
    L, M = peters_loads(a, b, U, ρ, bbar, θ, hdot, θdot, λ)
    # return portion of inputs that is not dependent on the state rates
    return SVector(θdot, 0, 0, L, M)
end

function get_input_state_jacobian(aero::PetersFiniteState, stru::TypicalSection,
    u, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # extract model constants
    bbar = aero.b
    # create zero row vector with length nλ
    zλ = zero(aero.b')
    # compute jacobian sub-matrices
    Jds = @SMatrix [0 0 0 1; 0 0 0 0; 0 0 0 0]
    Jda = vcat(zλ, zλ, zλ)
    Jrs = hcat(
        peters_loads_h(b, U, ρ),
        peters_loads_θ(b, U, ρ),
        peters_loads_hdot(),
        peters_loads_θdot(a, b, U, ρ)
        )
    Jra = peters_loads_λ(b, ρ, bbar)
    # return jacobian
    return vcat(hcat(Jds, Jda), hcat(Jrs, Jra))
end

# TODO: Parameter jacobian

# --- Internal Methods --- #
peters_lhs(a, b, Abar, cbar, dhdot, dθdot, dλ) = Abar*dλ
peters_rhs(a, b, U, cbar, θdot, hddot, θddot,  λ) = cbar*(hddot + U*θdot + (b/2-a*b)*θddot) - U/b*λ
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
