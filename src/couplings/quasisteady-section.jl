"""
    couple_models(aero::QuasiSteady, stru::TypicalSection)

Create an aerostructural model using a quasi-steady aerodynamics model and a
two-degree of freedom typical section model.  This model introduces the
freestream velocity ``U`` and air density ``\\rho`` as additional parameters.
"""
couple_models(aero::QuasiSteady, stru::TypicalSection)

# --- traits --- #
inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = 2

inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = 2

inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = 2

# --- methods --- #

function get_inputs(aero::QuasiSteady{0}, stru::TypicalSection, q, p, t)
    # extract state variables
    h, θ, hdot, θdot = q
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u = U
    v = U*θ
    # calculate aerodynamic loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # return inputs
    return SVector(L, M)
end

function get_inputs(aero::QuasiSteady{1}, stru::TypicalSection, q, p, t)
    # extract state variables
    h, θ, hdot, θdot = q
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # return inputs
    return SVector(L, M)
end

function get_inputs(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract state variables
    h, θ, hdot, θdot = q
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads
    L, M = quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, ω)
    # return inputs
    return SVector(L, M)
end

function get_input_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, ρ)
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady0_jacobian(a, b, ρ, a0, U)
end

function get_input_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady1_jacobian(a, b, ρ, a0, U)
end

function get_input_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady2_jacobian(a, b, ρ, a0, U)
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::QuasiSteady{0}, stru::TypicalSection,
    dq, q, p, t)

    return @SVector zeros(2)
end

function get_inputs_from_state_rates(aero::QuasiSteady{1}, stru::TypicalSection,
    dq, q, p, t)

    return @SVector zeros(2)
end

function get_inputs_from_state_rates(aero::QuasiSteady{2}, stru::TypicalSection,
    dq, q, p, t)
    # extract state rates
    dh, dθ, dhdot, dθdot = dq
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    # calculate aerodynamic loads
    L, M = quasisteady2_rate_loads(a, b, ρ, vdot, ωdot)
    # return inputs
    return SVector(L, M)
end
