"""
    Peters{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``d = \\begin{bmatrix}
u & \\dot{v} & \\dot{\\theta} & \\ddot{\\theta}\\end{bmatrix}^T`` and parameters
``p_a = \\begin{bmatrix} a & b & \\rho & a_0 & \\alpha_0 \\end{bmatrix}^T``
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
number_of_parameters(::Type{<:Peters}) = 5
inplaceness(::Type{<:Peters}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}) = Constant()
state_jacobian_type(::Type{<:Peters}) = Linear()
input_jacobian_type(::Type{<:Peters}) = Nonlinear()

# --- Methods --- #

function get_rates(model::Peters{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract inputs
    u, vdot, θdot, θddot = d
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return peters_rhs(a, b, cbar, u, vdot, θdot, θddot, λ)
end

get_mass_matrix(model::Peters) = model.A

# --- Performance Overloads --- #

function get_state_jacobian(model::Peters, λ, d, p, t)
    # extract inputs
    u, vdot, θdot, θddot = d
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return peters_state_jacobian(b, cbar, u)
end

function get_input_jacobian(model::Peters{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract inputs
    u, vdot, θdot, θddot = d
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return peters_input_jacobian(a, b, cbar, u, θdot, λ)
end

# --- Typical Section Coupling --- #

# --- traits --- #

inplaceness(::Type{<:Peters}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}, ::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{<:Peters}, ::Type{TypicalSection}) = 1

# --- methods --- #

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    s, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iλ = SVector{N}(1:N)
    iq = SVector{4}(N+1:N+4)
    # separate aerodynamic and structural states
    λ = s[iλ]
    q = s[iq]
    # extract structural state variables
    h, θ, hdot, θdot = q
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # calculate local vertical freestream velocity
    v = -u*θ - hdot
    # extract model constants
    bbar = aero.b
    # calculate aerodynamic loads
    L, M = peters_state_loads(a, b, ρ, a0, α0, bbar, u, v, θdot, λ)
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, 0, θdot, 0, L, M)
end

function get_input_mass_matrix(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 -1]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = hcat(
        zeros(SVector{2,TF}),
        zeros(SVector{2,TF}),
        -peters_loads_hddot(a, b, ρ),
        -peters_loads_θddot(a, b, ρ))
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # extract model constants
    bbar = aero.b
    # compute jacobian sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
    Jra = peters_loads_λ(a, b, ρ, a0, bbar, u)
    Jrs = hcat(
        peters_loads_h(),
        peters_loads_θ(a, b, ρ, a0, u),
        peters_loads_hdot(a, b, ρ, a0, u),
        peters_loads_θdot(a, b, ρ, a0, u)
        )
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    ds, s, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iλ = SVector{N}(1:N)
    iq = SVector{4}(N+1:N+4)
    # separate aerodynamic and structural state rates
    dλ = ds[iλ]
    dq = ds[iq]
    # extract structural state rates
    dh, dθ, dhdot, dθdot = dq
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    vdot = -dhdot
    θddot = dθdot
    # calculate aerodynamic loads
    L, M = peters_rate_loads(a, b, ρ, vdot, θddot)
    # return inputs
    return SVector(0, vdot, 0, θddot, L, M)
end

# --- Lifting Line Section Coupling --- #

# --- traits --- #

inplaceness(::Type{<:Peters}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Linear()
state_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{<:Peters}, ::Type{LiftingLineSection}) = 1

# --- methods --- #

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    s, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iλ = SVector{N}(1:N)
    iq = SVector{4}(N+1:N+6)
    # separate aerodynamic and structural states
    λ = s[iλ]
    q = s[iq]
    # extract structural state variables
    vx, vy, vz, ωx, ωy, ωz = q
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u = vx
    v = vz
    vdot = 0 # included in mass matrix
    θdot = ωy
    θddot = 0 # included in mass matrix
    # extract model constants
    bbar = aero.b
    # calculate aerodynamic loads
    L, M = peters_loads(a, b, ρ, a0, α0, bbar, u, v, vdot, θdot, θddot, λ)
    # forces and moments per unit span
    f = SVector(L, 0, 0)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return vcat(u, vdot, θdot, θddot, f, m)
end

function get_input_mass_matrix(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, s, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iλ = SVector{N}(1:N)
    iq = SVector{4}(N+1:N+6)
    # separate aerodynamic and structural states
    λ = s[iλ]
    q = s[iq]
    # extract structural state variables
    vx, vy, vz, ωx, ωy, ωz = q
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u = vx
    v = vz
    vdot = 0 # included in mass matrix
    θdot = ωy
    θddot = 0 # included in mass matrix
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 -1 0]
    Mra = zeros(SMatrix{2,N,TF})
    L_u, M_u = peters_loads_udot()
    L_v, M_v = peters_loads_vdot(a, b, ρ)
    L_θdot, M_θdot = peters_loads_θddot(a, b, ρ)
    Mrs = @SMatrix [-L_udot 0 -L_vdot 0 -L_θddot 0; 0 0 0 0 0 0; 0 0 0 0 0 0;
        0 0 0 0 0 0; -M_udot 0 -M_vdot 0 -M_θddot 0; 0 0 0 0 0 0]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, s, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iλ = SVector{N}(1:N)
    iq = SVector{4}(N+1:N+6)
    # separate aerodynamic and structural states
    λ = s[iλ]
    q = s[iq]
    # extract structural state variables
    vx, vy, vz, ωx, ωy, ωz = q
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u = vx
    v = vz
    vdot = 0 # included in mass matrix
    θdot = ωy
    θddot = 0 # included in mass matrix
    # extract model constants
    bbar = aero.b
    # compute jacobian sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0]
    Jra = peters_loads_λ(a, b, ρ, a0, bbar, u)
    L_u, M_u = peters_loads_u(a, b, ρ, a0, α0, bbar, u, v, vdot, θdot, θddot, λ)
    L_v, M_v = peters_loads_v(a, b, ρ, a0, u)
    L_θdot, M_θdot = peters_loads_θdot(a, b, ρ, a0, u)
    Jrs = @SMatrix [L_u 0 L_v 0 L_θdot 0; 0 0 0 0 0 0; 0 0 0 0 0 0;
        0 0 0 0 0 0; M_u 0 M_v 0 M_θdot 0; 0 0 0 0 0 0]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Peters, stru::LiftingLineSection,
    ds, s, p, t)
    # extract state rates
    dvx, dvy, dvz, dωx, dωy, dωz = ds
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    vdot = dvz
    θddot = dωy
    # calculate aerodynamic loads
    L, M = peters_rate_loads(a, b, ρ, vdot, θddot)
    # forces and moments per unit span
    f = SVector(L, 0, 0)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

# --- Internal Methods --- #

peters_lhs(Abar, dλ) = Abar*dλ

peters_rhs(a, b, cbar, u, vdot, θdot, θddot, λ) = cbar*(-vdot + u*θdot + (b/2-a*b)*θddot) - u/b*λ

peters_mass_matrix(Abar) = Abar

peters_state_jacobian(b, cbar, u) = -u/b*Diagonal(one.(cbar))

function peters_input_jacobian(a, b, cbar, u, θdot, λ)
    return hcat(cbar*θdot - 1/b*λ, -cbar, u*cbar, (b/2-a*b)*cbar)
end

function peters_loads(a, b, ρ, a0, α0, bbar, u, v, vdot, θdot, θddot, λ)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(-v + d*θdot - λ0 - u*α0) + tmp2*(-vdot/b + u/b*θdot - a*θddot)
    # moment at reference point
    M = -tmp2*(-vdot/2 + u*θdot + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function peters_state_loads(a, b, ρ, a0, α0, bbar, u, v, θdot, λ)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(-v + d*θdot - λ0 - u*α0) + tmp2*u/b*θdot
    # moment at reference point
    M = -tmp2*u*θdot + (b/2 + a*b)*L

    return SVector(L, M)
end

function peters_rate_loads(a, b, ρ, vdot, θddot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(-vdot/b - a*θddot)
    # moment at reference point
    M = -tmp*(-vdot/2 + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L

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

function peters_loads_u(a, b, ρ, a0, α0, bbar, u, v, vdot, θdot, θddot, λ)
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
    L_u = tmp1_u*(-v + d*θdot - λ0 - u*α0) - tmp1*α0 + tmp2/b*θdot
    # moment at reference point
    M_u = -tmp2*θdot + (b/2 + a*b)*L_u

    return SVector(L_u, M_u)
end

function peters_loads_v(a, b, ρ, a0, u)
    # lift at reference point
    L_v = -a0*ρ*u*b
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

peters_loads_udot() = SVector(0, 0)

function peters_loads_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_vdot = -tmp/b
    # moment at reference point
    M_vdot = tmp/2 + (b/2 + a*b)*L_vdot

    return SVector(L_vdot, M_vdot)
end

function peters_loads_θdot(a, b, ρ, a0, u)
    tmp = pi*ρ*b^3
    L_θdot = a0*ρ*u*b*(b/2 - a*b) + tmp*u/b
    M_θdot = -tmp*u + (b/2 + a*b)*L_θdot
    return SVector(L_θdot, M_θdot)
end

function peters_loads_θddot(a, b, ρ)
    tmp1 = pi*ρ*b^3
    tmp2 = b/2 + a*b
    L_θddot = -tmp1*a
    M_θddot = -tmp1*(b/8 - a*b/2) + tmp2*L_θddot
    return SVector(L_θddot, M_θddot)
end

peters_loads_h() = SVector(0, 0)

function peters_loads_θ(a, b, ρ, a0, u)
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end

function peters_loads_hdot(a, b, ρ, a0, u)
    L_hdot = a0*ρ*u*b
    M_hdot = (b/2 + a*b)*L_hdot
    return SVector(L_hdot, M_hdot)
end

function peters_loads_hddot(a, b, ρ)
    tmp1 = pi*ρ*b^3
    tmp2 = b/2 + a*b
    L_hddot = tmp1/b
    M_hddot = -tmp1/2 + tmp2*L_hddot
    return SVector(L_hddot, M_hddot)
end
