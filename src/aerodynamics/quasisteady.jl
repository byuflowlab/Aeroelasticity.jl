struct QuasiSteady{Order} <: NoStateModel end

# --- Constructors --- #

# zeroth order is steady state
Steady() = QuasiSteady{0}()

# default to using all terms
QuasiSteady() = QuasiSteady{2}()

# --- Model Traits --- #
number_of_parameters(::Type{QuasiSteady{Order}}) where Order = 4

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
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # calculate aerodynamic loads
    L, M = zero_order_loads(b, U, ρ, θ)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, neglecting apparent mass terms (except for moment calculation)
function get_inputs(aero::QuasiSteady{1}, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # calculate aerodynamic loads
    L, M = zero_order_loads(b, U, ρ, θ) .+
        first_order_loads(a, b, U, ρ, hdot, θdot)
    # return inputs
    return SVector(L, M)
end

# quasi-steady, including apparent mass terms
function get_inputs(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract state variables
    h, θ, hdot, θdot = u
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # calculate aerodynamic loads
    L, M = zero_order_loads(b, U, ρ, θ) .+
        first_order_loads(a, b, U, ρ, hdot, θdot) .+
        second_order_loads(a, b, U, ρ, θdot, 0.0, 0.0)
    # return inputs
    return SVector(L, M)
end

# steady
function get_input_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, u, p, t)
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return quasisteady0_jacobian(b, U, ρ)
end

# quasi-steady, neglecting acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, u, p, t)
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return quasisteady1_jacobian(a, b, U, ρ)
end

# quasi-steady, including acceleration terms
function get_input_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return quasisteady2_jacobian(a, b, U, ρ)
end

# quasi-steady, including acceleration terms
function get_input_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, u, p, t)
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, U, ρ)
end

# TODO: Parameter jacobian

# --- Internal --- #

function zero_order_loads(b, U, ρ, θ)
    L = 2*pi*ρ*b*U^2*θ
    M = 0 # quarter chord moment
    return SVector(L, M)
end

function first_order_loads(a, b, U, ρ, hdot, θdot)
    L = 2*pi*ρ*b*U*(hdot + (b/2 - a*b)*θdot)
    M = -pi*ρ*b^3*U*θdot # quarter chord moment
    return SVector(L, M)
end

function second_order_loads(a, b, U, ρ, θdot, hddot, θddot)
    L = pi*ρ*b^2*(hddot + U*θdot - b*a*θddot)
    M = -pi*ρ*b^3(hddot/2 + b*(1/8 - a/2)*θddot) # quarter chord moment
    return SVector(L, M)
end

function quasisteady0_jacobian(b, U, ρ)
    # calculate derivatives
    L_θ = 2*pi*ρ*b*U^2
    # return jacobian
    return @SMatrix [0 L_θ 0 0; 0 0 0 0]
end

function quasisteady1_jacobian(a, b, U, ρ)
    # calculate derivatives
    tmp = pi*ρ*b*U
    L_θ = 2*tmp*U
    L_hdot = 2*tmp
    L_θdot = 2*tmp*(b/2 - a*b)
    M_θdot = -tmp*b^2
    # return jacobian
    return @SMatrix [0 L_θ L_hdot L_θdot; 0 0 0 M_θdot]
end

function quasisteady2_jacobian(a, b, U, ρ)
    # calculate derivatives
    tmp = pi*ρ*b*U
    L_θ = 2*tmp*U
    L_hdot = 2*tmp
    L_θdot = 2*tmp*(b/2 - a*b) + tmp*b
    M_θdot = -tmp*b^2
    # return jacobian
    return @SMatrix [0 L_θ L_hdot L_θdot; 0 0 0 M_θdot]
end

function quasisteady2_mass_matrix(a, b, U, ρ)
    # calculate derivatives
    tmp = pi*ρ*b^2
    L_hddot = -pi*ρ*b^2
    L_θddot = pi*ρ*a*b^3
    M_hddot = pi*ρ*b^3*(1/2)
    M_θddot = pi*ρ*b^3*b*(1/8 - a/2)
    # return jacobian
    return @SMatrix [0 0 L_hddot L_θddot; 0 0 M_hddot M_θddot]
end
