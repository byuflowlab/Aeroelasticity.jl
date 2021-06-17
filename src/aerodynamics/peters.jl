"""
    Peters{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``d = \\begin{bmatrix}
u  & \\omega & \\dot{v} & \\dot{\\omega}\\end{bmatrix}^T`` and parameters
``p_a = \\begin{bmatrix} a & b & a_0 & \\alpha_0 \\end{bmatrix}^T``
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

# --- Typical Section Coupling --- #

"""
    couple_models(aero::Peters, stru::TypicalSection)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al. and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U`` and air density ``\\rho`` as additional
parameters.
"""
couple_models(aero::Peters, stru::TypicalSection)

# --- traits --- #

inplaceness(::Type{<:Peters}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}, ::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{<:Peters}, ::Type{TypicalSection}) = 2

# --- methods --- #

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    h, θ, hdot, θdot = s[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads
    L, M = peters_state_loads(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, 0, 0, L, M)
end

function get_input_mass_matrix(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    vdot_dhdot = 1
    ωdot_dθdot = 1
    # calculate aerodynamic loads
    L_dhdot, M_dhdot = peters_loads_vdot(a, b, ρ)
    L_dθdot, M_dθdot = peters_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 -vdot_dhdot 0; 0 0 0 -ωdot_dθdot]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = @SMatrix [0 0 -L_dhdot -L_dθdot; 0 0 -M_dhdot -M_dθdot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    ω_θdot = 1
    # calculate aerodynamic loads
    r_λ = peters_loads_λ(a, b, ρ, a0, bbar, U)
    L_θ, M_θ = peters_loads_θ(a, b, ρ, a0, U)
    L_hdot, M_hdot = peters_loads_v(a, b, ρ, a0, U)
    L_θdot, M_θdot = peters_loads_ω(a, b, ρ, a0, U)
    # construct sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 0; 0 0 0 ω_θdot; 0 0 0 0; 0 0 0 0]
    Jra = r_λ
    Jrs = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    ds, s, p, t) where {N,TF,SV,SA}
    # extract state rates
    dλ = ds[SVector{N}(1:N)]
    dh, dθ, dhdot, dθdot = ds[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream acceleration components
    vdot = dhdot
    ωdot = dθdot
    # calculate aerodynamic loads
    L, M = peters_rate_loads(a, b, ρ, vdot, ωdot)
    # return inputs
    return SVector(0, 0, vdot, ωdot, L, M)
end

# --- Lifting Line Section Coupling --- #

"""
    couple_models(aero::Peters, stru::LiftingLineSection)

Create an aerostructural model using a using the unsteady aerodynamic model
defined by Peters et al. and a lifting line section model.  The existence of this
coupling allows [`Peters`](@ref) to be used with [`LiftingLine`](@ref).  This
model introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::Peters, stru::LiftingLineSection)

# --- traits --- #

inplaceness(::Type{<:Peters}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Linear()
state_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{<:Peters}, ::Type{LiftingLineSection}) = 1

# --- methods --- #

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = s[SVector{6}(N+1:N+6)]
    # extract parameters
    a, b, a0, α0, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = peters_state_loads(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return vcat(u, ω, 0, 0, f, m)
end

function get_input_mass_matrix(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    vdot_dvz = 1
    ωdot_dωy = 1
    # calculate loads
    L_dvx, M_dvx = peters_loads_udot()
    L_dvz, M_dvz = peters_loads_vdot(a, b, ρ)
    L_dωy, M_dωy = peters_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0;
        0 0 -vdot_dvz 0 0 0; 0 0 0 0 -ωdot_dωy 0]
    Mra = zeros(SMatrix{6,N,TF})
    Mrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; -L_dvx 0 -L_dvz 0 -L_dωy 0;
        0 0 0 0 0 0; -M_dvx 0 -M_dvz 0 -M_dωy 0; 0 0 0 0 0 0]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = s[SVector{6}(N+1:N+6)]
    # extract parameters
    a, b, a0, α0, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    u_vx = 1
    ω_ωy = 1
    # compute loads
    out = peters_loads_λ(a, b, ρ, a0, bbar, u)
    L_λ, M_λ = out[1,:], out[2,:]
    L_vx, M_vx = peters_loads_u(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    L_vz, M_vz = peters_loads_v(a, b, ρ, a0, u)
    L_ωy, M_ωy = peters_loads_ω(a, b, ρ, a0, u)
    # construct submatrices
    Jda = zeros(SMatrix{4,N,TF})
    Jds = @SMatrix [u_vx 0 0 0 0 0; 0 0 0 0 ω_ωy 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
    Jra = vcat(zero(L_λ'), zero(L_λ'), L_λ', zero(M_λ'), M_λ', zero(M_λ'))
    Jrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_vx 0 L_vz 0 L_ωy 0;
        0 0 0 0 0 0; M_vx 0 M_vz 0 M_ωy 0; 0 0 0 0 0 0]
    # assemble jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, ds, s, p, t) where {N,TF,SV,SA}
    # extract state rates
    dλ = ds[SVector{N}(1:N)]
    dvx, dvy, dvz, dωx, dωy, dωz = ds[SVector{6}(N+1:N+6)]
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    vdot = dvz
    ωdot = dωy
    # calculate aerodynamic loads
    L, M = peters_rate_loads(a, b, ρ, vdot, ωdot)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(0, 0, vdot, ωdot, f, m)
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
