"""
    couple_models(aero::Peters, stru::TypicalSection, flap::LinearFlap)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al., a two-degree of freedom typical section model, and a linear
steady-state control surface model.  This model introduces the freestream
velocity ``U_\\infty``, air density ``\\rho_\\infty``, and flap deflection
``\\delta`` as additional parameters.
"""
couple_models(aero::Peters, stru::TypicalSection, flap::LinearFlap) = (aero, stru, flap)

# --- traits --- #

inplaceness(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{LinearFlap}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Linear()
state_jacobian_type(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{LinearFlap}) = Nonlinear()
number_of_parameters(::Type{<:Peters}, ::Type{TypicalSection}, ::Type{LinearFlap}) = 3

# --- methods --- #

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    flap::LinearFlap, s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    h, θ, hdot, θdot = s[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = U
    vbar = U*θ + hdot
    ω = θdot
    # calculate aerodynamic loads
    L, M = peters_state_loads(a, b, ρ, a0, α0, bbar, u, vbar, ω, λ)
    # add loads due to flap deflections
    L += ρ*U^2*b*clδ*δ
    M += 2*ρ*U^2*b^2*cmδ*δ
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, 0, 0, L, M)
end

function get_input_mass_matrix(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, flap::LinearFlap, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream velocity components
    vdot_dhdot = 1
    ωdot_dθdot = 1
    # calculate aerodynamic loads
    L_dhdot, M_dhdot = peters_loads_vdot(a, b, ρ)
    L_dθdot, M_dθdot = peters_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 -vdot_dhdot 0; 0 0 0 -ωdot_dθdot]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = @SMatrix [0 0 -L_dhdot -L_dθdot; 0 0 -M_dhdot -M_dθdot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::TypicalSection, flap::LinearFlap, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    ω_θdot = 1
    # calculate aerodynamic loads
    r_λ = peters_loads_λ(a, b, ρ, a0, bbar, U)
    L_θ, M_θ = peters_loads_θ(a, b, ρ, a0, U)
    L_hdot, M_hdot = peters_loads_v(a, b, ρ, a0, U)
    L_θdot, M_θdot = peters_loads_ω(a, b, ρ, a0, U)
    # construct sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 0; 0 0 0 ω_θdot; 0 0 0 0; 0 0 0 0]
    Jra = r_λ
    Jrs = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Peters{N,TF,SV,SA}, stru::TypicalSection,
    flap::LinearFlap, ds, s, p, t) where {N,TF,SV,SA}
    # extract state rates
    dλ = ds[SVector{N}(1:N)]
    dh, dθ, dhdot, dθdot = ds[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, clδ, cdδ, cmδ, U, ρ, δ = p
    # local freestream acceleration components
    vdot = dhdot
    ωdot = dθdot
    # calculate aerodynamic loads
    L, M = peters_rate_loads(a, b, ρ, vdot, ωdot)
    # return inputs
    return SVector(0, 0, vdot, ωdot, L, M)
end
