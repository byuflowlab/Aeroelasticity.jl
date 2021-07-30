"""
    couple_models(aero::Steady, stru::TypicalSection, flap::LinearFlap)

Create an aerostructural model using a steady aerodynamics model, a two-degree
of freedom typical section model, and a linear steady-state control surface
model.  This model introduces the freestream velocity ``U_\\infty``, air density
``\\rho_\\infty``, and flap deflection ``\\delta`` as additional parameters.
"""
couple_models(aero::Steady, stru::TypicalSection, flap::LinearFlap) = (aero, stru, flap)

"""
    couple_models(aero::QuasiSteady, stru::TypicalSection, flap::LinearFlap)

Create an aerostructural model using a quasi-steady aerodynamics model, a
two-degree of freedom typical section model, and a linear steady-state control
surface model.  This model introduces the freestream velocity ``U_\\infty``, air
density ``\\rho_\\infty``, and flap deflection ``\\delta`` as additional parameters.
"""
couple_models(aero::QuasiSteady, stru::TypicalSection, flap::LinearFlap) = (aero, stru, flap)

# --- traits --- #
inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = 3

inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Zeros()
state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = 3

inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = OutOfPlace()
mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Linear()
state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Nonlinear()
number_of_parameters(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{LinearFlap}) = 3

# --- methods --- #

function get_inputs(aero::QuasiSteady{0}, stru::TypicalSection, flap::LinearFlap,
    q, p, t)
    # extract state variables
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    u = U
    v = U*θ
    # calculate aerodynamic loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # add loads due to flap deflections
    L += ρ*U^2*b*clδ*δ
    M += 2*ρ*U^2*b^2*cmδ*δ
    # return inputs
    return SVector(L, M)
end

function get_inputs(aero::QuasiSteady{1}, stru::TypicalSection, flap::LinearFlap,
    q, p, t)
    # extract state variables
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # add loads due to flap deflections
    L += ρ*U^2*b*clδ*δ
    M += 2*ρ*U^2*b^2*cmδ*δ
    # return inputs
    return SVector(L, M)
end

function get_inputs(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract state variables
    h, θ, hdot, θdot = q
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads
    L, M = quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, ω)
    # add loads due to flap deflections
    L += ρ*U^2*b*clδ*δ
    M += 2*ρ*U^2*b^2*cmδ*δ
    # return inputs
    return SVector(L, M)
end

function get_input_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, ρ)
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # return jacobian
    return quasisteady0_jacobian(a, b, ρ, a0, U)
end

function get_input_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # return jacobian
    return quasisteady1_jacobian(a, b, ρ, a0, U)
end

function get_input_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
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
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    # calculate aerodynamic loads
    L, M = quasisteady2_rate_loads(a, b, ρ, vdot, ωdot)
    # return inputs
    return SVector(L, M)
end
