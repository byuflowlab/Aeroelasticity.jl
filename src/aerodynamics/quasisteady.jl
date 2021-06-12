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

# --- Model Traits --- #
number_of_parameters(::Type{<:QuasiSteady}) = 5

# --- Coupled Model Properties --- #

# steady assumption
inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = 1

# quasi-steady assumption, neglecting acceleration terms
inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = 1

# quasi-steady assumption, including acceleration terms
inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = 1

# --- Coupled Model Methods --- #

# steady
function get_inputs(aero::QuasiSteady{0}, stru::TypicalSection, s, p, t)
    # extract state variables
    h, θ, hdot, θdot = s
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    v = -u*θ
    # calculate aerodynamic loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, neglecting apparent mass terms (except for moment calculation)
function get_inputs(aero::QuasiSteady{1}, stru::TypicalSection, s, p, t)
    # extract state variables
    h, θ, hdot, θdot = s
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    v = -u*θ - hdot
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, θdot)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, including apparent mass terms
function get_inputs(aero::QuasiSteady{2}, stru::TypicalSection, s, p, t)
    # extract state variables
    h, θ, hdot, θdot = s
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # local vertical freestream velocity
    v = -u*θ - hdot
    vdot = 0 # defined with mass matrix
    θddot = 0 # defined with mass matrix
    # calculate aerodynamic loads
    L, M = quasisteady2_loads(a, b, ρ, a0, α0, u, v, vdot, θdot, θddot)
    # return inputs
    return SVector(L, M)
end

# steady
function get_input_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady0_jacobian(a, b, ρ, a0, u)
end

# quasi-steady, neglecting acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady1_jacobian(a, b, ρ, a0, u)
end

# quasi-steady, including acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady2_jacobian(a, b, ρ, a0, u)
end

# quasi-steady, including acceleration terms
function get_input_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, s, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, ρ, a0, α0, kh, kθ, m, Sθ, Iθ, u = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, ρ)
end

# TODO: Parameter jacobian

# --- Internal --- #

# steady state loads
function quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # lift at reference point
    L = -a0*ρ*b*u*v
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
    L = tmp1*(-v + d*θdot - u*α0) + tmp2*u/b*θdot
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
    L = tmp1*(-v + d*θdot - u*α0) + tmp2*(-vdot/b + u/b*θdot - a*θddot)
    # moment at reference point
    M = -tmp2*(-vdot/2 + u*θdot + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L

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
