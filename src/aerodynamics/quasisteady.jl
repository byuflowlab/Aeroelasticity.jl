"""
    QuasiSteady{Order} <: NoStateModel

2D quasi-steady aerodynamic model with parameters ``p_a = \\begin{bmatrix} a & b &
\\rho & a_0 & \alpha_0 \\end{bmatrix}^T``.
"""
struct QuasiSteady{Order} <: NoStateModel end

# --- Constructors --- #

"""
    Steady()

Initialize an object of type [`QuasiSteady`](@ref) which represents a 2D
steady aerodynamic model.
"""
Steady() = QuasiSteady{0}()

"""
    QuasiSteady()

Initialize an object of type [`QuasiSteady`](@ref) which represents a 2D
quasi-steady aerodynamic model.
"""
QuasiSteady() = QuasiSteady{2}()

# --- Traits --- #

number_of_parameters(::Type{<:QuasiSteady}) = 5

# --- Typical Section Coupling --- #

"""
    couple_models(aero::QuasiSteady, stru::TypicalSection)

Create an aerostructural model using a quasi-steady aerodynamics model and a
two-degree of freedom typical section model.  This model introduces the
freestream velocity ``U`` as an additional parameter.
"""
couple_models(aero::QuasiSteady, stru::TypicalSection)

# --- traits --- #
inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = 1

inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = 1

inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = 1

# --- methods --- #

function get_inputs(aero::QuasiSteady{0}, stru::TypicalSection, s, p, t)
    # extract state variables
    h, θ, hdot, θdot = s
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    v = u*θ
    # calculate aerodynamic loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # return inputs
    return SVector(L, M)
end

function get_inputs(aero::QuasiSteady{1}, stru::TypicalSection, s, p, t)
    # extract state variables
    h, θ, hdot, θdot = s
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    v = u*θ + hdot
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, θdot)
    # return inputs
    return SVector(L, M)
end

function get_inputs(aero::QuasiSteady{2}, stru::TypicalSection, s, p, t)
    # extract state variables
    h, θ, hdot, θdot = s
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    v = u*θ + hdot
    # calculate aerodynamic loads
    L, M = quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, θdot)
    # return inputs
    return SVector(L, M)
end

function get_input_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, ρ)
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady0_jacobian(a, b, ρ, a0, u)
end

function get_input_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady1_jacobian(a, b, ρ, a0, u)
end

function get_input_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady2_jacobian(a, b, ρ, a0, u)
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::QuasiSteady{0}, stru::TypicalSection,
    ds, s, p, t)

    return @SVector zeros(2)
end

function get_inputs_from_state_rates(aero::QuasiSteady{1}, stru::TypicalSection,
    ds, s, p, t)

    return @SVector zeros(2)
end

function get_inputs_from_state_rates(aero::QuasiSteady{2}, stru::TypicalSection,
    ds, s, p, t)
    # extract state rates
    dh, dθ, dhdot, dθdot = ds
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    vdot = dhdot
    θddot = dθdot
    # calculate aerodynamic loads
    L, M = quasisteady2_rate_loads(a, b, ρ, vdot, θddot)
    # return inputs
    return SVector(L, M)
end

# --- Lifting Line Section Coupling --- #

"""
    couple_models(aero::QuasiSteady, stru::LiftingLineSection)

Create an aerostructural model using a quasi-steady aerodynamics model and a
lifting line section model.  The existence of this coupling allows
[`QuasiSteady`](@ref) to be used with [`LiftingLine`](@ref).
"""
couple_models(aero::QuasiSteady, stru::LiftingLineSection)

# --- traits --- #

inplaceness(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = 1

inplaceness(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = 1

inplaceness(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Linear()
state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = 1

# --- methods --- #

function get_inputs(aero::QuasiSteady{0}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u, v = vx, vz
    # calculate loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # forces and moments per unit span
    f = SVector(L, 0, 0)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_inputs(aero::QuasiSteady{1}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u = vx
    v = vz
    θdot = ωy
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, θdot)
    # forces and moments per unit span
    f = SVector(L, 0, 0)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_inputs(aero::QuasiSteady{2}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u = vx
    v = vz
    θdot = ωy
    # calculate aerodynamic loads
    L, M = quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, θdot)
    # forces and moments per unit span
    f = SVector(L, 0, 0)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_input_mass_matrix(aero::QuasiSteady{2}, stru::LiftingLineSection,
    s, p, t)
    # extract parameters
    a, b, ρ, a0, α0 = p
    # calculate loads
    L_udot, M_udot = quasisteady2_udot()
    L_vdot, M_vdot = quasisteady2_vdot(a, b, ρ)
    L_θddot, M_θddot = quasisteady2_θddot(a, b, ρ)
    # return input mass matrix
    return @SMatrix [-L_udot 0 -L_vdot 0 -L_θddot 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;
        -M_udot 0 -M_vdot 0 -M_θddot 0; 0 0 0 0 0 0]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::QuasiSteady{0}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u, v = vx, vz
    # calculate loads
    L_u, M_u = quasisteady0_u(a, b, ρ, a0, v)
    L_v, M_v = quasisteady0_v(a, b, ρ, a0, u)
    # return inputs
    return @SMatrix [L_u 0 L_v 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;
        M_u 0 M_v 0 0 0; 0 0 0 0 0 0]
end

function get_input_state_jacobian(aero::QuasiSteady{1}, stru::LiftingLineSection,
    s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u = vx
    v = vz
    θdot = ωy
    # calculate loads
    L_u, M_u = quasisteady1_u(a, b, ρ, a0, α0, u, v, θdot)
    L_v, M_v = quasisteady1_v(a, b, ρ, a0, α0, u)
    L_θdot, M_θdot = quasisteady1_θdot(a, b, ρ, a0, u)
    # return inputs
    return @SMatrix [L_u 0 L_v 0 L_θdot 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;
        M_u 0 M_v 0 M_θdot 0; 0 0 0 0 0 0]
end

function get_input_state_jacobian(aero::QuasiSteady{2}, stru::LiftingLineSection,
    s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    u = vx
    v = vz
    vdot = 0 # included in mass matrix term
    θdot = ωy
    θddot = 0 # included in mass matrix term
    # calculate loads
    L_u, M_u = quasisteady2_u(a, b, ρ, a0, α0, u, v, θdot)
    L_v, M_v = quasisteady2_v(a, b, ρ, a0, α0, u)
    L_θdot, M_θdot = quasisteady2_θdot(a, b, ρ, a0, u)
    # return inputs
    return @SMatrix [L_u 0 L_v 0 L_θdot 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;
        M_u 0 M_v 0 M_θdot 0; 0 0 0 0 0 0]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::QuasiSteady{0}, stru::LiftingLineSection,
    ds, s, p, t)

    return @SVector zeros(6)
end

function get_inputs_from_state_rates(aero::QuasiSteady{1}, stru::LiftingLineSection,
    ds, s, p, t)
    
    return @SVector zeros(6)
end

function get_inputs_from_state_rates(aero::QuasiSteady{2}, stru::LiftingLineSection,
    ds, s, p, t)
    # extract state rates
    dvx, dvy, dvz, dωx, dωy, dωz = ds
    # extract parameters
    a, b, ρ, a0, α0 = p
    # extract relevant velocities
    vdot = dvz
    θddot = dωy
    # calculate aerodynamic loads
    L, M = quasisteady2_rate_loads(a, b, ρ, vdot, θddot)
    # forces and moments per unit span
    f = SVector(L, 0, 0)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

# --- Internal Methods --- #

# steady state loads
function quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # lift at reference point
    L = a0*ρ*b*u*v
    # moment at reference point
    M = (b/2 + a*b)*L

    return SVector(L, M)
end

# circulatory loads + added mass effects, neglecting accelerations
function quasisteady1_loads(a, b, ρ, a0, α0, u, v, θdot)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L = tmp1*(v + d*θdot - u*α0) + tmp2*u/b*θdot
    # moment at reference point
    M = -tmp2*u*θdot + (b/2 + a*b)*L

    return SVector(L, M)
end

# quasi-steady loads + added mass effects
function quasisteady2_loads(a, b, ρ, a0, α0, u, v, vdot, θdot, θddot)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L = tmp1*(v + d*θdot - u*α0) + tmp2*(vdot/b + u/b*θdot - a*θddot)
    # moment at reference point
    M = -tmp2*(vdot/2 + u*θdot + (b/8 - a*b/2)*θddot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, θdot)

    return quasisteady1_loads(a, b, ρ, a0, α0, u, v, θdot)
end

# quasi-steady loads from acceleration terms
function quasisteady2_rate_loads(a, b, ρ, vdot, θddot)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(vdot/b - a*θddot)
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function quasisteady0_jacobian(a, b, ρ, a0, u)
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    return @SMatrix [0 L_θ 0 0; 0 M_θ 0 0]
end

function quasisteady1_jacobian(a, b, ρ, a0, u)
    tmp1 = a0*ρ*u*b
    tmp2 = pi*ρ*b^3
    d1 = b/2 - a*b
    d2 = b/2 + a*b
    L_θ = tmp1*u
    L_hdot = tmp1
    L_θdot = tmp1*d1 + tmp2*u/b
    M_θ = d2*L_θ
    M_hdot = d2*L_hdot
    M_θdot = -tmp2*u + d2*L_θdot
    return @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
end

quasisteady2_jacobian(a, b, ρ, a0, u) = quasisteady1_jacobian(a, b, ρ, a0, u)

function quasisteady2_mass_matrix(a, b, ρ)
    # calculate derivatives
    tmp1 = pi*ρ*b^3
    tmp2 = b/2 + a*b
    L_hddot = -tmp1/b
    L_θddot = tmp1*a
    M_hddot = tmp1/2 + tmp2*L_hddot
    M_θddot = tmp1*(b/8 - a*b/2) + tmp2*L_θddot
    # return jacobian
    return @SMatrix [0 0 L_hddot L_θddot; 0 0 M_hddot M_θddot]
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

function quasisteady1_u(a, b, ρ, a0, α0, u, v, θdot)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    tmp1_u = a0*ρ*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L_u = tmp1_u*(v + d*θdot - u*α0) - tmp1*α0 + tmp2/b*θdot
    # moment at reference point
    M_u = -tmp2*θdot + (b/2 + a*b)*L_u

    return SVector(L_u, M_u)
end

function quasisteady1_v(a, b, ρ, a0, α0, u)
    # lift at reference point
    L_v = a0*ρ*u*b
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

function quasisteady1_θdot(a, b, ρ, a0, u)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L_θdot = tmp1*d + tmp2*u/b
    # moment at reference point
    M_θdot = -tmp2*u + (b/2 + a*b)*L_θdot

    return SVector(L_θdot, M_θdot)
end

quasisteady2_u(a, b, ρ, a0, α0, u, v, θdot) = quasisteady1_u(a, b, ρ, a0, α0, u, v, θdot)

quasisteady2_v(a, b, ρ, a0, α0, u) = quasisteady1_v(a, b, ρ, a0, α0, u)

quasisteady2_udot() = SVector(0, 0)

function quasisteady2_vdot(a, b, ρ)
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_vdot = tmp/b
    # moment at reference point
    M_vdot = -tmp/2 + (b/2 + a*b)*L_vdot

    return SVector(L_vdot, M_vdot)
end

quasisteady2_θdot(a, b, ρ, a0, u) = quasisteady1_θdot(a, b, ρ, a0, u)

function quasisteady2_θddot(a, b, ρ)
    tmp = pi*ρ*b^3
    L_θddot = -tmp*a
    M_θddot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_θddot
    return SVector(L_θddot, M_θddot)
end
