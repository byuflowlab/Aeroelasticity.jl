"""
    couple_models(aero::Peters, stru::LiftingLineSection)

Create an aerostructural model using a using the unsteady aerodynamic model
defined by Peters et al. and a lifting line section model.  The existence of this
coupling allows [`Peters`](@ref) to be used with [`LiftingLine`](@ref).  This
model introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::Peters, stru::LiftingLineSection) = (aero, stru)

# --- Traits --- #

number_of_additional_parameters(::Type{<:Peters}, ::Type{LiftingLineSection}) = 1

coupling_inplaceness(::Type{<:Peters}, ::Type{LiftingLineSection}) = OutOfPlace()

coupling_rate_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Linear()
coupling_state_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Nonlinear()
coupling_time_gradient_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Zeros()

# --- Methods --- #

function get_coupling_inputs(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    dx, x, p, t) where {N,TF,SV,SA}

    # extract rate variables
    dλ = dx[SVector{N}(1:N)]
    dvx, dvy, dvz, dωx, dωy, dωz = dx[SVector{6}(N+1:N+6)]

    # extract state variables
    λ = x[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = x[SVector{6}(N+1:N+6)]

    # extract parameters
    a, b, a0, α0, cd0, cm0, ρ = p

    # extract model constants
    bbar = aero.b

    # freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_accelerations(dvx, dvz, dωy)

    # aerodynamic loads
    Na, Aa, Ma = peters_loads(a, b, ρ, a0, α0, cd0, cm0, bbar, u, v, ω, vdot, ωdot, λ)

    # loads per unit span
    f = SVector(Aa, 0, Na)
    m = SVector(0, Ma, 0)

    # return inputs
    return SVector(u, ω, vdot, ωdot, f..., m...)
end

# --- Performance Overloads --- #

#TODO: State rate, state, and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Peters, stru::LiftingLineSection; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::Peters, stru::LiftingLineSection, padd)

    return (rho = padd[1],)
end
