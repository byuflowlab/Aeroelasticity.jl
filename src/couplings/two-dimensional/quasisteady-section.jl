"""
    couple_models(aero::Steady, stru::TypicalSection)

Create an aerostructural model using a steady aerodynamics model and a
two-degree of freedom typical section model.  This model introduces the
freestream velocity ``U`` and air density ``\\rho`` as additional parameters.
"""
couple_models(aero::Steady, stru::TypicalSection) = (aero, stru)

"""
    couple_models(aero::QuasiSteady, stru::TypicalSection)

Create an aerostructural model using a quasi-steady aerodynamics model and a
two-degree of freedom typical section model.  This model introduces the
freestream velocity ``U`` and air density ``\\rho`` as additional parameters.
"""
couple_models(aero::QuasiSteady, stru::TypicalSection) = (aero, stru)

# --- traits --- #

number_of_additional_parameters(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = OutOfPlace()
coupling_mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}) = Nonlinear()

number_of_additional_parameters(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = OutOfPlace()
coupling_mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}) = Nonlinear()

number_of_additional_parameters(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = OutOfPlace()
coupling_mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Linear()
coupling_state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}) = Nonlinear()

# --- methods --- #

function get_coupling_inputs(aero::QuasiSteady{0}, stru::TypicalSection, q, p, t)
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

function get_coupling_inputs(aero::QuasiSteady{1}, stru::TypicalSection, q, p, t)
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

function get_coupling_inputs(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
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

function get_coupling_mass_matrix(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady2_mass_matrix(a, b, ρ)
end

# --- performance overloads --- #

function get_coupling_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady0_jacobian(a, b, ρ, a0, U)
end

function get_coupling_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady1_jacobian(a, b, ρ, a0, U)
end

function get_coupling_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection, q, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # return jacobian
    return quasisteady2_jacobian(a, b, ρ, a0, U)
end

# --- unit testing methods --- #

function get_coupling_inputs_using_state_rates(aero::QuasiSteady{0}, stru::TypicalSection,
    dq, q, p, t)

    return @SVector zeros(2)
end

function get_coupling_inputs_using_state_rates(aero::QuasiSteady{1}, stru::TypicalSection,
    dq, q, p, t)

    return @SVector zeros(2)
end

function get_coupling_inputs_using_state_rates(aero::QuasiSteady{2}, stru::TypicalSection,
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

# --- convenience methods --- #

function set_additional_parameters!(padd, aero::QuasiSteady, stru::TypicalSection; U, rho)

    padd[1] = U
    padd[2] = rho

    return padd
end

function separate_additional_parameters(aero::QuasiSteady, stru::TypicalSection, padd)

    return (U = padd[1], rho = padd[2])
end

# --- Plotting --- #

@recipe function f(aero::QuasiSteady, stru::TypicalSection, x, y, p, t)

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    h, θ, hdot, θdot = x

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end
