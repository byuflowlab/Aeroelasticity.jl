struct QuasiSteady{Order} <: NoStateModel end

# --- Constructors --- #

# zeroth order is steady state
Steady() = QuasiSteady{0}()

# default to using all terms
QuasiSteady() = QuasiSteady{2}()

# --- Model Traits --- #
number_of_parameters(::Type{<:QuasiSteady}) = 4

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
    a, b, U, ρ = p
    # calculate aerodynamic loads
    L, M = zero_order_loads(a, b, U, ρ, θ)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, neglecting apparent mass terms (except for moment calculation)
function get_inputs(aero::QuasiSteady{1}, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # calculate aerodynamic loads
    L, M = zero_order_loads(a, b, U, ρ, θ) .+
        first_order_loads(a, b, U, ρ, hdot, θdot)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, including apparent mass terms
function get_inputs(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # calculate aerodynamic loads
    L, M = zero_order_loads(a, b, U, ρ, θ) .+
        first_order_loads(a, b, U, ρ, hdot, θdot) .+
        second_order_loads(a, b, U, ρ, θdot, 0.0, 0.0)
    # return inputs
    return SVector(L, M)
end

# steady
function get_input_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # return jacobian
    return quasisteady0_jacobian(a, b, U, ρ)
end

# quasi-steady, neglecting acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # return jacobian
    return quasisteady1_jacobian(a, b, U, ρ)
end

# quasi-steady, including acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # return jacobian
    return quasisteady2_jacobian(a, b, U, ρ)
end

# quasi-steady, including acceleration terms
function get_input_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract aerodynamic parameters
    a, b, U, ρ = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, U, ρ)
end

# TODO: Parameter jacobian

# --- Internal --- #

function zero_order_loads(a, b, U, ρ, θ)
    L = 2*pi*ρ*b*U^2*θ
    M = (b/2 + a*b)*L
    return SVector(L, M)
end

function first_order_loads(a, b, U, ρ, hdot, θdot)
    L = 2*pi*ρ*b*U*(hdot + (b/2 - a*b)*θdot)
    M = -pi*ρ*b^3*U*θdot + (b/2 + a*b)*L
    return SVector(L, M)
end

function second_order_loads(a, b, U, ρ, θdot, hddot, θddot)
    L = pi*ρ*b^2*(hddot + U*θdot - b*a*θddot)
    M = -pi*ρ*b^3*(hddot/2 + b*(1/8 - a/2)*θddot) + (b/2 + a*b)*L
    return SVector(L, M)
end

function quasisteady0_jacobian(a, b, U, ρ)
    L_θ = 2*pi*ρ*b*U^2
    M_θ = (b/2 + a*b)*L_θ
    return @SMatrix [0 L_θ 0 0; 0 M_θ 0 0]
end

function quasisteady1_jacobian(a, b, U, ρ)
    # calculate derivatives
    tmp1 = pi*ρ*b*U
    tmp2 = b/2 + a*b
    L_θ = 2*tmp1*U
    L_hdot = 2*tmp1
    L_θdot = 2*tmp1*(b/2 - a*b)
    M_θ = tmp2*L_θ
    M_hdot = tmp2*L_hdot
    M_θdot = -tmp1*b^2 + tmp2*L_θdot
    # return jacobian
    return @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
end

function quasisteady2_jacobian(a, b, U, ρ)
    # calculate derivatives
    tmp1 = pi*ρ*b*U
    tmp2 = b/2 + a*b
    L_θ = 2*tmp1*U
    L_hdot = 2*tmp1
    L_θdot = 2*tmp1*(b/2 - a*b) + tmp1*b
    M_θ = tmp2*L_θ
    M_hdot = tmp2*L_hdot
    M_θdot = -tmp1*b^2 + tmp2*L_θdot
    # return jacobian
    return @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
end

function quasisteady2_mass_matrix(a, b, U, ρ)
    # calculate derivatives
    tmp1 = pi*ρ*b^2
    tmp2 = b/2 + a*b
    L_hddot = -tmp1
    L_θddot = tmp1*a*b
    M_hddot = tmp1*b/2 + tmp2*L_hddot
    M_θddot = tmp1*(b^2/8 - a*b^2/2) + tmp2*L_θddot
    # return jacobian
    return @SMatrix [0 0 L_hddot L_θddot; 0 0 M_hddot M_θddot]
end
