"""
    couple_models(aero::Wagner, stru::TypicalSection, flap::SimpleFlap)

Create an aerostructural model using an unsteady aerodynamic model based on
Wagner's function, a two-degree of freedom typical section model, and a linear
steady-state control surface model.  This model introduces the freestream
velocity ``U_\\infty``, air density ``\\rho_\\infty``, and flap deflection
``\\delta`` as additional parameters.
"""
couple_models(aero::Wagner, stru::TypicalSection, flap::SimpleFlap) = (aero, stru, flap)

# --- Traits --- #

number_of_additional_parameters(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = 3

coupling_inplaceness(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = OutOfPlace()

coupling_rate_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Linear()
coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()
coupling_time_gradient_type(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Zeros()

# --- Methods --- #

function get_coupling_inputs(aero::Wagner, stru::TypicalSection, flap::SimpleFlap, dx, x, p, t)
    # extract rate variables
    dλ1, dλ2, dh, dθ, dhdot, dθdot = dx
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate aerodynamic loads (except contribution from state rates)
    La, Ma = wagner_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # loads due to flap deflections
    Lf, Df, Mf = simpleflap_loads(b, u, ρ, clδ, cdδ, cmδ, δ)
    # loads per unit span
    L = La + Lf
    M = Ma + Mf
    # return inputs
    return SVector(u, v, ω, L, M)
end

# --- Performance Overloads --- #

function get_coupling_rate_jacobian(aero::Wagner, stru::TypicalSection,
    flap::SimpleFlap, dx, x, p, t)
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # calculate loads
    L_hddot, M_hddot = wagner_loads_vdot(a, b, ρ)
    L_θddot, M_θddot = wagner_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0]
    Mrs = @SMatrix [0 0 L_hddot L_θddot; 0 0 M_hddot M_θddot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

function get_coupling_state_jacobian(aero::Wagner, stru::TypicalSection,
    flap::SimpleFlap, dx, x, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u_θ, v_θ, ω_θ = section_velocities_θ(U)
    u_hdot, v_hdot, ω_hdot = section_velocities_hdot()
    u_θdot, v_θdot, ω_θdot = section_velocities_θdot()
    # calculate loads
    r_λ = wagner_loads_λ(a, b, ρ, a0, U)
    L_θ, M_θ = wagner_loads_θ(a, b, ρ, a0, C1, C2, U)
    L_hdot, M_hdot = wagner_loads_v(a, b, ρ, a0, C1, C2, U)
    L_θdot, M_θdot = wagner_loads_ω(a, b, ρ, a0, C1, C2, U)
    # compute jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [0 0 0 0; 0 v_θ v_hdot 0; 0 0 0 ω_θdot]
    Jra = r_λ
    Jrs = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# TODO: Parameter jacobian

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Wagner, stru::TypicalSection,
    flap::SimpleFlap; U, rho, delta)

    padd[1] = U
    padd[2] = rho
    padd[3] = delta

    return padd
end

function separate_additional_parameters(aero::Wagner, stru::TypicalSection,
    flap::SimpleFlap, padd)

    return (U = padd[1], rho = padd[2], delta = padd[3])
end

# --- Plotting --- #

@recipe function f(aero::Wagner, stru::TypicalSection, flap::SimpleFlap, dx, x, y, p, t)

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end
