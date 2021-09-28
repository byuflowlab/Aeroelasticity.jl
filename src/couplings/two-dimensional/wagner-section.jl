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
coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Linear()
coupling_parameter_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Nonlinear()
coupling_time_gradient_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Zeros()

# --- Methods --- #

function get_coupling_inputs(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract rate variables
    dλ1, dλ2, dh, dθ, dhdot, dθdot = dx
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    N, A, M = wagner_loads(a, b, ρ, a0, α0, cd0, cm0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # lift is approximately normal force
    L = N
    # return inputs
    return SVector(u, v, ω, L, M)
end

# --- Performance Overloads --- #

function get_coupling_rate_jacobian(aero::Wagner, stru::TypicalSection, p)
    # extract parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # calculate loads
    N_hddot, A_hddot, M_hddot = wagner_loads_vdot(a, b, ρ)
    N_θddot, A_θddot, M_θddot = wagner_loads_ωdot(a, b, ρ)
    # construct submatrices
    Jyaa = @SMatrix [0 0; 0 0; 0 0]
    Jyas = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Jysa = @SMatrix [0 0; 0 0]
    Jyss = @SMatrix [0 0 N_hddot N_θddot; 0 0 M_hddot M_θddot]
    # assemble rate jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_state_jacobian(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    u_θ, v_θ, ω_θ = section_velocities_θ(U)
    u_hdot, v_hdot, ω_hdot = section_velocities_hdot()
    u_θdot, v_θdot, ω_θdot = section_velocities_θdot()
    # calculate loads
    r_λ = wagner_loads_λ(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    N_λ, A_λ, M_λ = r_λ[1,:], r_λ[2,:], r_λ[3,:]
    N_θ, M_θ = wagner_loads_θ(a, b, ρ, a0, C1, C2, U)
    N_hdot, A_hdot, M_hdot = wagner_loads_v(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    N_θdot, A_θdot, M_θdot = wagner_loads_ω(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # compute jacobian sub-matrices
    Jyaa = @SMatrix [0 0; 0 0; 0 0]
    Jyas = @SMatrix [0 u_θ u_hdot u_θdot; 0 v_θ v_hdot v_θdot; 0 ω_θ ω_hdot ω_θdot]
    Jysa = vcat(N_λ', M_λ')
    Jyss = @SMatrix [0 N_θ N_hdot N_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_parameter_jacobian(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract rate variables
    dλ1, dλ2, dh, dθ, dhdot, dθdot = dx
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    u_U, v_U, ω_U = section_velocities_U(θ)
    # calculate loads
    N_a, A_a, M_a = wagner_loads_a(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    N_b, A_b, M_b = wagner_loads_b(a, b, ρ, a0, α0, cd0, cm0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    N_a0, A_a0, M_a0 = wagner_loads_a0(a, b, ρ, α0, C1, C2, u, v, ω, λ1, λ2)
    N_α0, A_α0, M_α0 = wagner_loads_α0(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    N_cd0, A_cd0, M_cd0 = wagner_loads_cd0(b, ρ, u)
    N_cm0, A_cm0, M_cm0 = wagner_loads_cm0(b, ρ, u)

    N_u, A_u, M_u = wagner_loads_u(a, b, ρ, a0, α0, cd0, cm0, C1, C2, u, v, ω, λ1, λ2)
    N_v, A_v, M_v = wagner_loads_v(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    N_U = N_u * u_U + N_v * v_U
    M_U = M_u * u_U + M_v * v_U

    N_ρ, A_ρ, M_ρ = wagner_loads_ρ(a, b, a0, α0, cd0, cm0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
        # compute jacobian sub-matrices
    Jyaa = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
    Jyas = @SMatrix [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
    Jyap = @SMatrix [u_U 0; v_U 0; ω_U 0]
    Jysa = @SMatrix [N_a N_b N_a0 N_α0 N_cd0 N_cm0; M_a M_b M_a0 M_α0 M_cd0 M_cm0]
    Jyss = @SMatrix [0 0 0 0 0; 0 0 0 0 0]
    Jysp = @SMatrix [N_U N_ρ; M_U M_ρ]
    # return jacobian
    return [Jyaa Jyas Jyap; Jysa Jyss Jysp]
end

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Wagner, stru::TypicalSection; U, rho)

    padd[1] = U
    padd[2] = rho

    return padd
end

function separate_additional_parameters(aero::Wagner, stru::TypicalSection, padd)

    return (U = padd[1], rho = padd[2])
end

# --- Plotting --- #

@recipe function f(aero::Wagner, stru::TypicalSection, dx, x, y, p, t)

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x

    # extract parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end

# --- Internal Methods --- #

function wagner_loads_θ(a, b, ρ, a0, C1, C2, U)
    ϕ0 = 1 - C1 - C2
    N_θ = a0*ρ*b*U^2*ϕ0
    M_θ = (b/2 + a*b)*N_θ
    return SVector(N_θ, M_θ)
end
