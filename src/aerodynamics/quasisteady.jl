struct QuasiSteady{Order} <: NoStateModel end

# --- Constructors --- #

# zeroth order is steady state
Steady() = QuasiSteady{0}()

# default to using all terms
QuasiSteady() = QuasiSteady{2}()

# --- Model Traits --- #
number_of_parameters(::Type{<:QuasiSteady}) = 6

# --- Coupled Model Properties --- #

# steady assumption
inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Varying()

# quasi-steady assumption, neglecting acceleration terms
inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Varying()

# quasi-steady assumption, including acceleration terms
inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Varying()
state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Varying()

# --- Coupled Model Methods --- #

# steady
function get_inputs(aero::QuasiSteady{0}, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # calculate aerodynamic loads
    L, M = quasisteady0_loads(a, b, U, ρ, a0, α0, θ)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, neglecting apparent mass terms (except for moment calculation)
function get_inputs(aero::QuasiSteady{1}, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, U, ρ, a0, α0, θ, hdot, θdot)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, including apparent mass terms
function get_inputs(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # calculate aerodynamic loads
    L, M = quasisteady2_loads(a, b, U, ρ, a0, α0, θ, hdot, θdot, 0.0, 0.0)
    # return inputs
    return SVector(L, M)
end

# steady
function get_input_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # return jacobian
    return quasisteady0_jacobian(a, b, U, ρ, a0)
end

# quasi-steady, neglecting acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # return jacobian
    return quasisteady1_jacobian(a, b, U, ρ, a0)
end

# quasi-steady, including acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # return jacobian
    return quasisteady2_jacobian(a, b, U, ρ, a0)
end

# quasi-steady, including acceleration terms
function get_input_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ, a0, α0 = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, U, ρ)
end

# TODO: Parameter jacobian

# --- Internal --- #

# steady state loads
function quasisteady0_loads(a, b, U, ρ, a0, α0, θ)
    # circulatory load factor
    tmp = a0*ρ*U^2*b
    # effective angle of attack (assuming small angles)
    α = θ
    # lift at reference point
    L = tmp*θ
    # moment at reference point
    M = (b/2 + a*b)*L

    return SVector(L, M)
end

# circulatory loads + added mass effects, neglecting accelerations
function quasisteady1_loads(a, b, U, ρ, a0, α0, θ, hdot, θdot)
    # circulatory load factor
    tmp1 = a0*ρ*U*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L = tmp1*(U*θ + hdot + d*θdot - U*α0) + tmp2*U/b*θdot
    # moment at reference point
    M = -tmp2*U*θdot + (b/2 + a*b)*L

    return SVector(L, M)
end

# quasi-steady loads + added mass effects
function quasisteady2_loads(a, b, U, ρ, a0, α0, θ, hdot, θdot, hddot, θddot)
    # circulatory load factor
    tmp1 = a0*ρ*U*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # lift at reference point
    L = tmp1*(U*θ + hdot + d*θdot - U*α0) + tmp2*(hddot/b + U/b*θdot - a*θddot)
    # moment at reference point
    M = -tmp2*(hddot/2 + U*θdot + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L

    return SVector(L, M)
end

function quasisteady0_jacobian(a, b, U, ρ, a0)
    L_θ = a0*ρ*U^2*b
    M_θ = (b/2 + a*b)*L_θ
    return @SMatrix [0 L_θ 0 0; 0 M_θ 0 0]
end

function quasisteady1_jacobian(a, b, U, ρ, a0)
    tmp1 = a0*ρ*U*b
    tmp2 = pi*ρ*b^3
    d1 = b/2 - a*b
    d2 = b/2 + a*b
    L_θ = tmp1*U
    L_hdot = tmp1
    L_θdot = tmp1*d1 + tmp2*U/b
    M_θ = d2*L_θ
    M_hdot = d2*L_hdot
    M_θdot = -tmp2*U + d2*L_θdot
    return @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
end

quasisteady2_jacobian(a, b, U, ρ, a0) = quasisteady1_jacobian(a, b, U, ρ, a0)

function quasisteady2_mass_matrix(a, b, U, ρ)
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
