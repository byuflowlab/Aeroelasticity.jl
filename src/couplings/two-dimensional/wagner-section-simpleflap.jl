"""
    couple_models(aero::Wagner, stru::TypicalSection, flap::SimpleFlap)

Create an aerostructural model using an unsteady aerodynamic model based on
Wagner's function, a two-degree of freedom typical section model, and a linear
steady-state control surface model.  This model introduces the freestream
velocity ``U_\\infty``, air density ``\\rho_\\infty``, and flap deflection
``\\delta`` as additional parameters.
"""
couple_models(aero::Wagner, stru::TypicalSection, flap::SimpleFlap) = (aero, stru, flap)

# --- traits --- #

number_of_additional_parameters(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = 3
coupling_inplaceness(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = OutOfPlace()
coupling_mass_matrix_type(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Linear()
coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}, ::Type{SimpleFlap}) = Nonlinear()

# --- methods --- #

function get_coupling_inputs(aero::Wagner, stru::TypicalSection, flap::SimpleFlap, s, p, t)
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = s
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads (except contribution from state rates)
    L, M = wagner_state_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # add loads due to flap deflections
    L += ρ*U^2*b*clδ*δ
    M += 2*ρ*U^2*b^2*cmδ*δ
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, v, ω, L, M)
end

function get_coupling_mass_matrix(aero::Wagner, stru::TypicalSection,
    flap::SimpleFlap, s, p, t)
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # calculate loads
    L_hddot, M_hddot = wagner_loads_vdot(a, b, ρ)
    L_θddot, M_θddot = wagner_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0]
    Mrs = @SMatrix [0 0 -L_hddot -L_θddot; 0 0 -M_hddot -M_θddot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_coupling_state_jacobian(aero::Wagner, stru::TypicalSection,
    flap::SimpleFlap, u, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    v_θ = U
    v_hdot = 1
    ω_θdot = 1
    # calculate loads
    r_λ = wagner_loads_λ(a, b, ρ, a0, U)
    L_h, M_h = wagner_loads_h()
    L_θ, M_θ = wagner_loads_θ(a, b, ρ, a0, C1, C2, U)
    L_hdot, M_hdot = wagner_loads_v(a, b, ρ, a0, C1, C2, U)
    L_θdot, M_θdot = wagner_loads_ω(a, b, ρ, a0, C1, C2, U)
    # compute jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [0 0 0 0; 0 v_θ v_hdot 0; 0 0 0 ω_θdot]
    Jra = r_λ
    Jrs = @SMatrix [L_h L_θ L_hdot L_θdot; M_h M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_coupling_inputs_using_state_rates(aero::Wagner, stru::TypicalSection,
    flap::SimpleFlap, ds, s, p, t)
    # extract state rates
    dλ1, dλ2, dh, dθ, dhdot, dθdot = ds
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    # calculate aerodynamic loads
    L, M = wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # return inputs
    return SVector(0, 0, 0, L, M)
end

# --- convenience methods --- #

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

@recipe function f(aero::Wagner, stru::TypicalSection, flap::SimpleFlap, x, y, p, t)

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
