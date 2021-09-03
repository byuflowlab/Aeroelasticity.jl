"""
    couple_models(aero::Peters, stru::TypicalSection)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al. and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::Peters, stru::TypicalSection) = (aero, stru)

# --- traits --- #

number_of_additional_parameters(::Type{<:Peters}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{<:Peters}, ::Type{TypicalSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Linear()
coupling_state_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Nonlinear()

# --- methods --- #

function get_coupling_inputs(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    dx, x, p, t) where {N,TF,SV,SA}
    # extract rate variables
    dλ = dx[SVector{N}(1:N)]
    dh, dθ, dhdot, dθdot = dx[SVector{4}(N+1:N+4)]
    # extract state variables
    λ = x[SVector{N}(1:N)]
    h, θ, hdot, θdot = x[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    L, M = peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    # return inputs
    return SVector(u, ω, vdot, ωdot, L, M)
end

# --- performance overloads --- #

function get_coupling_rate_jacobian(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, dx, x, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    udot_dhdot, vdot_dhdot, ωdot_dhdot = section_accelerations_dhdot()
    udot_dθdot, vdot_dθdot, ωdot_dθdot = section_accelerations_dθdot()
    # calculate aerodynamic loads
    L_dhdot, M_dhdot = peters_loads_vdot(a, b, ρ)
    L_dθdot, M_dθdot = peters_loads_ωdot(a, b, ρ)
    # construct submatrices
    Jyaa = zeros(SMatrix{4,N,TF})
    Jyas = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 vdot_dhdot vdot_dθdot; 0 0 ωdot_dhdot ωdot_dθdot]
    Jysa = zeros(SMatrix{2,N,TF})
    Jyss = @SMatrix [0 0 L_dhdot L_dθdot; 0 0 M_dhdot M_dθdot]
    # return jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, x, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u_θ, v_θ, ω_θ = section_velocities_θ(U)
    u_hdot, v_hdot, ω_hdot = section_velocities_hdot()
    u_θdot, v_θdot, ω_θdot = section_velocities_θdot()
    # calculate aerodynamic loads
    r_λ = peters_loads_λ(a, b, ρ, a0, bbar, U)
    L_θ, M_θ = peters_loads_θ(a, b, ρ, a0, U)
    L_hdot, M_hdot = peters_loads_v(a, b, ρ, a0, U)
    L_θdot, M_θdot = peters_loads_ω(a, b, ρ, a0, U)
    # construct sub-matrices
    Jyaa = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jyas = @SMatrix [0 0 0 0; 0 ω_θ ω_hdot ω_θdot; 0 0 0 0; 0 0 0 0]
    Jysa = r_λ
    Jyss = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_parameter_jacobian(aero::Wagner, stru::TypicalSection, dx, x, p, t)
    # extract rate variables
    dλ = dx[SVector{N}(1:N)]
    dh, dθ, dhdot, dθdot = dx[SVector{4}(N+1:N+4)]
    # extract state variables
    λ = x[SVector{N}(1:N)]
    h, θ, hdot, θdot = x[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u_U, v_U, ω_U = section_velocities_U(U)
    # calculate loads
    L_a, M_a = peters_loads_a(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    L_b, M_b = peters_loads_b(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    L_a0, M_a0 = peters_loads_a0(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    L_α0, M_α0 = peters_loads_α0(a, b, ρ, a0, u)
    L_U, M_U = peters_loads_u(a, b, ρ, a0, α0, C1, C2, u, v, ωdot, λ1, λ2)
    L_ρ, M_ρ = peters_loads_ρ(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
        # compute jacobian sub-matrices
    Jyaa = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Jyas = @SMatrix [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
    Jyap = @SMatrix [u_U 0; ω_U 0; 0 0; 0 0]
    Jysa = @SMatrix [L_a L_b L_a0 L_α0; M_a M_b M_a0 M_α0]
    Jyss = @SMatrix [0 0 0 0 0; 0 0 0 0 0]
    Jysp = @SMatrix [L_U L_ρ; M_U M_ρ]
    # return jacobian
    return [Jyaa Jyas Jyap; Jysa Jyss Jysp]
end

# --- convenience methods --- #

function set_additional_parameters!(padd, aero::Peters, stru::TypicalSection; U, rho)

    padd[1] = U
    padd[2] = rho

    return padd
end

function separate_additional_parameters(aero::Peters, stru::TypicalSection, padd)

    return (U = padd[1], rho = padd[2])
end

# --- plotting --- #

@recipe function f(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    x, y, p, t) where {N,TF,SV,SA}

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    λ = x[SVector{N}(1:N)]
    h, θ, hdot, θdot = x[SVector{4}(N+1:N+4)]

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end

# --- Internal Methods --- #

function peters_loads_θ(a, b, ρ, a0, u)
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end
