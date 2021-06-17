"""
    Wagner{TF} <: AbstractModel

Aerodynamic model based on Wagner's function with state variables ``\\lambda =
\\begin{bmatrix} \\lambda_1 & \\lambda_2 \\end{bmatrix}^T``, inputs ``d =
\\begin{bmatrix} u & v & \\omega \\end{bmatrix}^T``
and parameters ``p_a = \\begin{bmatrix} a & b & a_0 & \alpha_0 \\end{bmatrix}^T``
"""
struct Wagner{TF} <: AbstractModel
    C1::TF
    C2::TF
    eps1::TF
    eps2::TF
end

# --- Constructors --- #

"""
    Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

Initialize an object of type [`Wagner`](@ref)
"""
function Wagner(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)
    return Wagner(C1, C2, eps1, eps2)
end

# --- Traits --- #

number_of_states(::Type{<:Wagner}) = 2
number_of_inputs(::Type{<:Wagner}) = 3
number_of_parameters(::Type{<:Wagner}) = 4
inplaceness(::Type{<:Wagner}) = OutOfPlace()
mass_matrix_type(::Type{<:Wagner}) = Identity()
state_jacobian_type(::Type{<:Wagner}) = Linear()
input_jacobian_type(::Type{<:Wagner}) = Nonlinear()

# --- Methods --- #

function get_rates(model::Wagner, λ, d, p, t) where {N,TF,SV,SA}
    # extract states
    λ1, λ2 = λ
    # extract inputs
    u, v, ω = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # calculate rates
    return wagner_rates(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)
end

# --- Performance Overloads --- #

function get_state_jacobian(model::Wagner, λ, d, p, t)
    # extract inputs
    u, v, ω = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    ε1 = model.eps1
    ε2 = model.eps2
    # jacobian with respect to aerodynamic states
    return wagner_state_jacobian(a, b, ε1, ε2, u)
end

function get_input_jacobian(model::Wagner, λ, d, p, t)
    # extract states
    λ1, λ2 = λ
    # extract inputs
    u, v, ω = d
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2
    # return jacobian
    return wagner_input_jacobian(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)
end

# --- Unit Testing Methods --- #

get_lhs(::Wagner, dλ, λ, d, p, t) = dλ

# --- Typical Section Coupling --- #

"""
    couple_models(aero::Wagner, stru::TypicalSection)

Create an aerostructural model using an unsteady aerodynamic model based on
Wagner's function and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U`` and air density ``\\rho`` as additional
parameters.
"""
couple_models(aero::Wagner, stru::TypicalSection)

# --- traits --- #

inplaceness(::Type{<:Wagner}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{<:Wagner}, ::Type{TypicalSection}) = 2

# --- methods --- #

function get_inputs(aero::Wagner, stru::TypicalSection, s, p, t)
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = s
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads (except contribution from state rates)
    L, M = wagner_state_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, v, ω, L, M)
end

function get_input_mass_matrix(aero::Wagner, stru::TypicalSection, s, p, t)
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # calculate loads
    L_hddot, M_hddot = wagner_loads_vdot(a, b, ρ)
    L_θddot, M_θddot = wagner_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0]
    Mrs = @SMatrix [0 0 -L_hddot -L_θddot; 0 0 -M_hddot -M_θddot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Wagner, stru::TypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    v_θ = U
    v_hdot = 1
    ω_θdot = 1
    # calculate loads
    r_λ = wagner_loads_λ(a, b, ρ, a0, U)
    L_h, M_h = wagner_loads_h()
    L_θ, M_θ = wagner_loads_θ(a, b, ρ, a0, C1, C2, U)
    L_hdot, M_hdot = wagner_loads_v(a, b, ρ, a0, C1, C2, U)
    L_θdot, M_θdot = wagner_loads_ω(a, b, ρ, a0, C1, C2, U)
    # compute jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [0 0 0 0; 0 v_θ v_hdot 0; 0 0 0 ω_θdot]
    Jra = r_λ
    Jrs = @SMatrix [L_h L_θ L_hdot L_θdot; M_h M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Wagner, stru::TypicalSection,
    ds, s, p, t)
    # extract state rates
    dλ1, dλ2, dh, dθ, dhdot, dθdot = ds
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    # calculate aerodynamic loads
    L, M = wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # return inputs
    return SVector(0, 0, 0, L, M)
end

# --- Lifting Line Section Coupling --- #

"""
    couple_models(aero::Wagner, stru::LiftingLineSection)

Create an aerostructural model using an unsteady aerodynamic model based
on Wagner's function and a lifting line section model.  The existence of this
coupling allows [`Wagner`](@ref) to be used with [`LiftingLine`](@ref).  This
model introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::Wagner, stru::LiftingLineSection)

# --- traits --- #

inplaceness(::Type{<:Wagner}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Linear()
state_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{<:Wagner}, ::Type{LiftingLineSection}) = 1

# --- methods --- #

function get_inputs(aero::Wagner, stru::LiftingLineSection, s, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # calculate aerodynamic loads
    L, M = wagner_state_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return vcat(u, v, ω, f, m)
end

function get_input_mass_matrix(aero::Wagner, stru::LiftingLineSection, s, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # calculate loads
    L_dvx, M_dvx = wagner_loads_udot()
    L_dvz, M_dvz = wagner_loads_vdot(a, b, ρ)
    L_dωy, M_dωy = wagner_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0; 0 0; 0 0; 0 0; 0 0]
    Mrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; -L_dvx 0 -L_dvz 0 -L_dωy 0;
        0 0 0 0 0 0; -M_dvx 0 -M_dvz 0 -M_dωy 0; 0 0 0 0 0 0]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Wagner, stru::LiftingLineSection, s, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    u_vx = 1
    v_vz = 1
    ω_ωy = 1
    # model constants
    C1 = aero.C1
    C2 = aero.C2
    # calculate loads
    out = wagner_loads_λ(a, b, ρ, a0, u)
    L_λ, M_λ = out[1,:], out[2,:]
    L_vx, M_vx = wagner_loads_u(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    L_vz, M_vz = wagner_loads_v(a, b, ρ, a0, C1, C2, u)
    L_ωy, M_ωy = wagner_loads_ω(a, b, ρ, a0, C1, C2, u)
    # compute input jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [u_vx 0 0 0 0 0; 0 0 v_vz 0 0 0; 0 0 0 0 ω_ωy 0]
    Jra = vcat(zero(L_λ'), zero(L_λ'), L_λ', zero(M_λ'), M_λ', zero(M_λ'))
    Jrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_vx 0 L_vz 0 L_ωy 0;
        0 0 0 0 0 0; M_vx 0 M_vz 0 M_ωy 0; 0 0 0 0 0 0]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Wagner, stru::LiftingLineSection,
    ds, s, p, t)
    # extract state rates
    dλ1, dλ2, dvx, dvy, dvz, dωx, dωy, dωz = ds
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    vdot = dvz
    ωdot = dωy
    # calculate aerodynamic loads
    L, M = wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(0, 0, 0, f, m)
end

# --- Internal Methods --- #

function wagner_rates(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)
    λ1dot = -ε1*u/b*λ1 + C1*ε1*u/b*(v + (1/2-a)*b*ω - u*α0)
    λ2dot = -ε2*u/b*λ2 + C2*ε2*u/b*(v + (1/2-a)*b*ω - u*α0)
    return SVector(λ1dot, λ2dot)
end

wagner_state_jacobian(a, b, ε1, ε2, u) = @SMatrix [-ε1*u/b 0; 0 -ε2*u/b]

function wagner_input_jacobian(a, b, α0, C1, C2, ε1, ε2, u, v, ω, λ1, λ2)

    tmp1 = v/b + (1/2-a)*ω - 2*α0*u/b
    λ1dot_u = -ε1/b*λ1 + C1*ε1*tmp1
    λ2dot_u = -ε2/b*λ2 + C2*ε2*tmp1

    tmp2 = u/b
    λ1dot_v = C1*ε1*tmp2
    λ2dot_v = C2*ε2*tmp2

    tmp3 = u*(1/2-a)
    λ1dot_ω = C1*ε1*tmp3
    λ2dot_ω = C2*ε2*tmp3

    return @SMatrix [λ1dot_u λ1dot_v λ1dot_ω; λ2dot_u λ2dot_v λ2dot_ω]
end

function wagner_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*(vdot/b + u/b*ω - a*ωdot)
    # moment at reference point
    M = -tmp2*(vdot/2 + u*ω + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function wagner_state_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L

    return SVector(L, M)
end

function wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(vdot/b - a*ωdot)
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function wagner_loads_λ(a, b, ρ, a0, u)
    tmp1 = a0*ρ*u*b
    tmp2 = (b/2 + a*b)*tmp1
    return @SMatrix [tmp1 tmp1; tmp2 tmp2]
end
wagner_loads_λdot() = @SMatrix [0 0; 0 0]

function wagner_loads_u(a, b, ρ, a0, α0, C1, C2, u, v, ωdot, λ1, λ2)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    tmp1_u = a0*ρ*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_u = tmp1_u*((v + d*ωdot - u*α0)*ϕ0 + λ1 + λ2) - tmp1*α0*ϕ0 + tmp2/b*ωdot
    # moment at reference point
    M_u = -tmp2*ωdot + (b/2 + a*b)*L_u

    return SVector(L_u, M_u)
end

function wagner_loads_v(a, b, ρ, a0, C1, C2, u)
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_v = a0*ρ*u*b*ϕ0
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

function wagner_loads_ω(a, b, ρ, a0, C1, C2, u)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_ω = tmp1*d*ϕ0 + tmp2*u/b
    # moment at reference point
    M_ω = -tmp2*u + (b/2 + a*b)*L_ω

    return SVector(L_ω, M_ω)
end

wagner_loads_udot() = SVector(0, 0)

function wagner_loads_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_vdot = tmp/b
    # moment at reference point
    M_vdot = -tmp/2 + (b/2 + a*b)*L_vdot

    return SVector(L_vdot, M_vdot)
end

function wagner_loads_ωdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_ωdot = -a*tmp
    # moment at reference point
    M_ωdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_ωdot

    return SVector(L_ωdot, M_ωdot)
end

wagner_loads_h() = SVector(0, 0)

function wagner_loads_θ(a, b, ρ, a0, C1, C2, u)
    ϕ0 = 1 - C1 - C2
    L_θ = a0*ρ*b*u^2*ϕ0
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end
