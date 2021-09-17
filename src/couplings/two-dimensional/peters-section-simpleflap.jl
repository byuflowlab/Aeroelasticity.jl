"""
    couple_models(aero::Peters, stru::TypicalSection, flap::SimpleFlap)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al., a two-degree of freedom typical section model, and a linear
steady-state control surface model.  This model introduces the freestream
velocity ``U_\\infty``, air density ``\\rho_\\infty``, and flap deflection
``\\delta`` as additional parameters.
"""
couple_models(aero::Peters, stru::TypicalSection, flap::SimpleFlap) = (aero, stru, flap)

# --- Traits --- #

number_of_additional_parameters(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = 3

coupling_inplaceness(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = OutOfPlace()

coupling_rate_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Linear()
coupling_state_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_time_gradient_type(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Zeros()

# --- Methods --- #

function get_coupling_inputs(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    flap::SimpleFlap, dx, x, p, t) where {N,TF,SV,SA}

    # extract rate variables
    dλ = dx[SVector{N}(1:N)]
    dh, dθ, dhdot, dθdot = dx[SVector{4}(N+1:N+4)]

    # extract state variables
    λ = x[SVector{N}(1:N)]
    h, θ, hdot, θdot = x[SVector{4}(N+1:N+4)]

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p

    # extract model constants
    bbar = aero.b

    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)

    # aerodynamic loads
    La, Ma = peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    # loads due to flap deflections
    Lf, Df, Mf = simpleflap_loads(b, u, ρ, clδ, cdδ, cmδ, δ)

    # loads per unit span
    L = La + Lf
    M = Ma + Mf

    # return inputs
    return SVector(u, ω, vdot, ωdot, L, M)
end

# --- Performance Overloads --- #

# TODO: State rate, state, and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Peters, stru::TypicalSection,
    flap::SimpleFlap; U, rho, delta)

    padd[1] = U
    padd[2] = rho
    padd[3] = delta

    return padd
end

function separate_additional_parameters(aero::Peters, stru::TypicalSection,
    flap::SimpleFlap, padd)

    return (U = padd[1], rho = padd[2], delta = padd[3])
end

# --- Plotting --- #

@recipe function f(aero::Peters{N,TF,SV,SA}, stru::TypicalSection, flap::SimpleFlap,
    dx, x, y, p, t) where {N,TF,SV,SA}

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    λ = x[SVector{N}(1:N)]
    h, θ, hdot, θdot = x[SVector{4}(N+1:N+4)]

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end
