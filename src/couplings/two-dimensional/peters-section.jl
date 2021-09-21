"""
    couple_models(aero::Peters, stru::TypicalSection)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al. and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::Peters, stru::TypicalSection) = (aero, stru)

# --- Traits --- #

number_of_additional_parameters(::Type{<:Peters}, ::Type{TypicalSection}) = 2

coupling_inplaceness(::Type{<:Peters}, ::Type{TypicalSection}) = OutOfPlace()

coupling_rate_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Linear()
coupling_state_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}) = Nonlinear()
coupling_time_gradient_type(::Type{<:Peters}, ::Type{TypicalSection}) = Zeros()

# --- Methods --- #

function get_coupling_inputs(aero::Peters{Nλ,TF,SV,SA}, stru::TypicalSection,
    dx, x, p, t) where {Nλ,TF,SV,SA}
    # extract rate variables
    dλ = dx[SVector{Nλ}(1:Nλ)]
    dh, dθ, dhdot, dθdot = dx[SVector{4}(Nλ+1:Nλ+4)]
    # extract state variables
    λ = x[SVector{Nλ}(1:Nλ)]
    h, θ, hdot, θdot = x[SVector{4}(Nλ+1:Nλ+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    N, A, M = peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return SVector(u, ω, vdot, ωdot, L, M)
end

# --- Performance Overloads --- #

function get_coupling_rate_jacobian(aero::Peters{Nλ,TF,SV,SA},
    stru::TypicalSection, dx, x, p, t) where {Nλ,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    udot_dhdot, vdot_dhdot, ωdot_dhdot = section_accelerations_dhdot()
    udot_dθdot, vdot_dθdot, ωdot_dθdot = section_accelerations_dθdot()
    # calculate aerodynamic loads
    N_vdot, A_vdot, M_vdot = peters_loads_vdot(a, b, ρ)
    N_ωdot, A_ωdot, M_ωdot = peters_loads_ωdot(a, b, ρ)
    # propagate derivatives
    N_dhdot = N_vdot * vdot_dhdot + N_ωdot * ωdot_dhdot
    N_dθdot = N_vdot * vdot_dθdot + N_ωdot * ωdot_dθdot
    M_dhdot = M_vdot * vdot_dhdot + M_ωdot * ωdot_dhdot
    M_dθdot = M_vdot * vdot_dθdot + M_ωdot * ωdot_dθdot
    # construct submatrices
    Jyaa = zeros(SMatrix{4,Nλ,TF})
    Jyas = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 vdot_dhdot vdot_dθdot; 0 0 ωdot_dhdot ωdot_dθdot]
    Jysa = zeros(SMatrix{2,Nλ,TF})
    Jyss = @SMatrix [0 0 N_dhdot N_dθdot; 0 0 M_dhdot M_dθdot]
    # return jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_state_jacobian(aero::Peters{Nλ,TF,SV,SA},
    stru::TypicalSection, dx, x, p, t) where {Nλ,TF,SV,SA}
    # extract state variables
    λ = x[SVector{Nλ}(1:Nλ)]
    h, θ, hdot, θdot = x[SVector{4}(Nλ+1:Nλ+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    u_θ, v_θ, ω_θ = section_velocities_θ(U)
    u_hdot, v_hdot, ω_hdot = section_velocities_hdot()
    u_θdot, v_θdot, ω_θdot = section_velocities_θdot()
    # calculate aerodynamic loads
    r_λ = peters_loads_λ(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    N_λ, A_λ, M_λ = r_λ[1,:], r_λ[2,:], r_λ[3,:]
    N_θ, M_θ = peters_loads_θ(a, b, ρ, a0, U)
    N_hdot, A_hdot, M_hdot = peters_loads_v(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    N_θdot, A_θdot, M_θdot = peters_loads_ω(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    # construct sub-matrices
    Jyaa = zeros(SMatrix{4,Nλ,TF}) # d(d)/d(dλ)
    Jyas = @SMatrix [0 0 0 0; 0 ω_θ ω_hdot ω_θdot; 0 0 0 0; 0 0 0 0]
    Jysa = vcat(N_λ', M_λ')
    Jyss = @SMatrix [0 N_θ N_hdot N_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jyaa Jyas; Jysa Jyss]
end

function get_coupling_parameter_jacobian(aero::Peters{Nλ,TF,SV,SA},
    stru::TypicalSection, dx, x, p, t) where {Nλ,TF,SV,SA}
    # extract rate variables
    dλ = dx[SVector{Nλ}(1:Nλ)]
    dh, dθ, dhdot, dθdot = dx[SVector{4}(Nλ+1:Nλ+4)]
    # extract state variables
    λ = x[SVector{Nλ}(1:Nλ)]
    h, θ, hdot, θdot = x[SVector{4}(Nλ+1:Nλ+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    u_U, v_U, ω_U = section_velocities_U(θ)
    # calculate loads
    N_a, A_a, M_a = peters_loads_a(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    N_b, A_b, M_b = peters_loads_b(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
    N_a0, A_a0, M_a0 = peters_loads_a0(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    N_α0, A_α0, M_α0 = peters_loads_α0(a, b, ρ, a0, α0, bbar, u, v, ω, λ)

    N_u, A_u, M_u = peters_loads_u(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    N_v, A_v, M_v = peters_loads_v(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    N_U = N_u * u_U + N_v * v_U
    M_U = M_u * u_U + M_v * v_U

    N_ρ, A_ρ, M_ρ = peters_loads_ρ(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)
        # compute jacobian sub-matrices
    Jyaa = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    Jyas = @SMatrix [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
    Jyap = @SMatrix [u_U 0; ω_U 0; 0 0; 0 0]
    Jysa = @SMatrix [N_a N_b N_a0 N_α0; M_a M_b M_a0 M_α0]
    Jyss = @SMatrix [0 0 0 0 0; 0 0 0 0 0]
    Jysp = @SMatrix [N_U N_ρ; M_U M_ρ]
    # return jacobian
    return [Jyaa Jyas Jyap; Jysa Jyss Jysp]
end

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Peters, stru::TypicalSection; U, rho)

    padd[1] = U
    padd[2] = rho

    return padd
end

function separate_additional_parameters(aero::Peters, stru::TypicalSection, padd)

    return (U = padd[1], rho = padd[2])
end

# --- Plotting --- #

@recipe function f(aero::Peters{Nλ,TF,SV,SA}, stru::TypicalSection,
    dx, x, y, p, t) where {Nλ,TF,SV,SA}

    framestyle --> :origin
    grid --> false
    xlims --> (-1.0, 1.0)
    ylims --> (-1.5, 1.5)
    label --> @sprintf("t = %6.3f", t)

    # extract state variables
    λ = x[SVector{Nλ}(1:Nλ)]
    h, θ, hdot, θdot = x[SVector{4}(Nλ+1:Nλ+4)]

    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p

    xplot = [-(0.5 + a*b)*cos(θ),    (0.5 - a*b)*cos(θ)]
    yplot = [ (0.5 + a*b)*sin(θ)-h, -(0.5 - a*b)*sin(θ)-h]

    return xplot, yplot
end

# --- Internal Methods --- #

function peters_loads_θ(a, b, ρ, a0, u)
    N_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*N_θ
    return SVector(N_θ, M_θ)
end
