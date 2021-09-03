"""
    couple_models(aero::Wagner, stru::TypicalSection)

Create an aerostructural model using an unsteady aerodynamic model based on
Wagner's function and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::Wagner, stru::TypicalSection) = (aero, stru)

# --- Traits --- #

number_of_additional_parameters(::Type{<:Wagner}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{<:Wagner}, ::Type{TypicalSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Constant()
coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Nonlinear()

# --- Methods --- #

function get_coupling_inputs(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract rate variables
    dλ1, dλ2, dh, dθ, dhdot, dθdot = dx
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    L, M = wagner_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # return inputs
    return SVector(u, v, ω, L, M)
end

# --- Performance Overloads --- #

function get_coupling_rate_jacobian(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # calculate loads
    L_hddot, M_hddot = wagner_loads_vdot(a, b, ρ)
    L_θddot, M_θddot = wagner_loads_ωdot(a, b, ρ)
    # construct submatrices
    Jyaa = @SMatrix [0 0; 0 0; 0 0]
    Jyas = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Jysa = @SMatrix [0 0; 0 0]
    Jyss = @SMatrix [0 0 L_hddot L_θddot; 0 0 M_hddot M_θddot]
    # assemble rate jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_state_jacobian(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
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
    Jyaa = @SMatrix [0 0; 0 0; 0 0]
    Jyas = @SMatrix [0 u_θ u_hdot u_θdot; 0 v_θ v_hdot v_θdot; 0 ω_θ ω_hdot ω_θdot]
    Jysa = r_λ
    Jyss = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_parameter_jacobian(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract rate variables
    dλ1, dλ2, dh, dθ, dhdot, dθdot = dx
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u_U, v_U, ω_U = section_velocities_U(U)
    # calculate loads
    L_a, M_a = wagner_loads_a(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    L_b, M_b = wagner_loads_b(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    L_a0, M_a0 = wagner_loads_a0(a, b, ρ, α0, C1, C2, u, v, ω, λ1, λ2)
    L_α0, M_a0 = wagner_loads_α0(a, b, ρ, a0, C1, C2, u)
    L_U, M_U = wagner_loads_u(a, b, ρ, a0, α0, C1, C2, u, v, ωdot, λ1, λ2)
    L_ρ, M_ρ = wagner_loads_ρ(a, b, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
        # compute jacobian sub-matrices
    Jyaa = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Jyas = @SMatrix [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
    Jyap = @SMatrix [u_U 0; v_U 0; ω_U 0]
    Jysa = @SMatrix [L_a L_b L_a0 L_α0; M_a M_b M_a0 M_α0]
    Jyss = @SMatrix [0 0 0 0 0; 0 0 0 0 0]
    Jysp = @SMatrix [L_U L_ρ; M_U M_ρ]
    # return jacobian
    return [Jyaa Jyas Jyap; Jysa Jyss Jysp]
end

# --- convenience methods --- #

function set_additional_parameters!(padd, aero::Wagner, stru::TypicalSection; U, rho)

    padd[1] = U
    padd[2] = rho

    return padd
end

function separate_additional_parameters(aero::Wagner, stru::TypicalSection, padd)

    return (U = padd[1], rho = padd[2])
end

# --- Plotting --- #

@recipe function f(aero::Wagner, stru::TypicalSection, x, y, p, t)

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end

# --- Internal Methods --- #

function wagner_loads_θ(a, b, ρ, a0, C1, C2, U)
    ϕ0 = 1 - C1 - C2
    L_θ = a0*ρ*b*U^2*ϕ0
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end
