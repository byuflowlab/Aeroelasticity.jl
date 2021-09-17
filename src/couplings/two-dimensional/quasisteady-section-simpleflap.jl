"""
    couple_models(aero::Steady, stru::TypicalSection, flap::SimpleFlap)

Create an aerostructural model using a steady aerodynamics model, a
two-degree of freedom typical section model, and a linear steady-state control
surface model.  This model introduces the freestream velocity ``U``, air density
``\\rho``, and control surface deflection ``\\delta`` as additional parameters.
"""
function couple_models(aero::Steady, stru::TypicalSection, flap::SimpleFlap)
    return (aero, stru, flap)
end

"""
    couple_models(aero::QuasiSteady, stru::TypicalSection)

Create an aerostructural model using a quasi-steady aerodynamics model, a
two-degree of freedom typical section model, and a linear steady-state control
surface model.  This model introduces the freestream velocity ``U``, air density
``\\rho``, and control surface deflection ``\\delta`` as additional parameters.
"""
function couple_models(aero::QuasiSteady, stru::TypicalSection, flap::SimpleFlap)
    return (aero, stru, flap)
end

# --- traits --- #

number_of_additional_parameters(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = 3
coupling_inplaceness(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_time_gradient_type(::Type{QuasiSteady{0}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Zeros()

number_of_additional_parameters(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = 3
coupling_inplaceness(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_time_gradient_type(::Type{QuasiSteady{1}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Zeros()

number_of_additional_parameters(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = 3
coupling_inplaceness(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Linear()
coupling_state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_time_gradient_type(::Type{QuasiSteady{2}}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Zeros()

# --- methods --- #

function get_coupling_inputs(aero::QuasiSteady{0}, stru::TypicalSection, flap::SimpleFlap,
    dx, x, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    u, v = section_steady_velocities(U, θ)
    # calculate aerodynamic loads
    La, Ma = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # loads due to flap deflections
    Lf, Df, Mf = simpleflap_loads(b, u, ρ, clδ, cdδ, cmδ, δ)
    # total loads
    L = La + Lf
    M = Ma + Mf
    # return inputs
    return SVector(L, M)
end

function get_coupling_inputs(aero::QuasiSteady{1}, stru::TypicalSection, flap::SimpleFlap,
    dx, x, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    # calculate aerodynamic loads
    La, Ma = quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # loads due to flap deflections
    Lf, Df, Mf = simpleflap_loads(b, u, ρ, clδ, cdδ, cmδ, δ)
    # total loads
    L = La + Lf
    M = Ma + Mf
    # return inputs
    return SVector(L, M)
end

function get_coupling_inputs(aero::QuasiSteady{2}, stru::TypicalSection, flap::SimpleFlap,
    dx, x, p, t)
    # extract rate variables
    dh, dθ, dhdot, dθdot = dx
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # aerodynamic loads
    La, Ma = quasisteady2_loads(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)
    # loads due to flap deflections
    Lf, Df, Mf = simpleflap_loads(b, u, ρ, clδ, cdδ, cmδ, δ)
    # total loads
    L = La + Lf
    M = Ma + Mf
    # return inputs
    return SVector(L, M)
end

# --- Performance Overloads --- #

function get_coupling_state_jacobian(aero::QuasiSteady{0}, stru::TypicalSection,
    flap::SimpleFlap, dx, x, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # return jacobian
    return quasisteady0_section_state_jacobian(a, b, ρ, a0, U)
end

function get_coupling_state_jacobian(aero::QuasiSteady{1}, stru::TypicalSection,
    flap::SimpleFlap, dx, x, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # return jacobian
    return quasisteady1_section_state_jacobian(a, b, ρ, a0, U)
end

function get_coupling_state_jacobian(aero::QuasiSteady{2}, stru::TypicalSection,
    flap::SimpleFlap, dx, x, p, t)
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # return jacobian
    return quasisteady2_section_state_jacobian(a, b, ρ, a0, U)
end

# TODO: State rate and parameter jacobians

# --- convenience methods --- #

function set_additional_parameters!(padd, aero::QuasiSteady, stru::TypicalSection,
    flap::SimpleFlap; U, rho, delta)

    padd[1] = U
    padd[2] = rho
    padd[3] = delta

    return padd
end

function separate_additional_parameters(aero::QuasiSteady, stru::TypicalSection,
    flap::SimpleFlap, padd)

    return (U = padd[1], rho = padd[2], delta = padd[3])
end

# --- plotting --- #

@recipe function f(aero::QuasiSteady, stru::TypicalSection, flap::SimpleFlap,
    dx, x, y, p, t)

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    h, θ, hdot, θdot = x

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end
