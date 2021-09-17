"""
    couple_models(aero::Peters, stru::LiftingLineSection, flap::SimpleFlap,
        ctrl::LiftingLineSectionControl)

Create an aerostructural model using a using the unsteady aerodynamic model
defined by Peters et al, a lifting line aerodynamic model, and a linear steady-state
control surface model.  The existence of this coupling allows [`Peters`](@ref)
and [`SimpleFlap`](@ref) to be used with [`LiftingLine`](@ref) and
[`LiftingLineFlaps`](@ref).  This model introduces the freestream air density
``\\rho`` as an additional parameter.
"""
function couple_models(aero::Peters, stru::LiftingLineSection, flap::SimpleFlap,
    ctrl::LiftingLineSectionControl)

    return (aero, stru, flap, ctrl)
end

# --- Traits --- #

function number_of_additional_parameters(::Type{<:Peters}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return 1
end

function coupling_inplaceness(::Type{<:Peters}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function coupling_rate_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Linear()
end

function coupling_state_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_parameter_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_time_gradient_type(::Type{<:Peters}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

# --- Methods --- #

function get_coupling_inputs(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, dx, x, p, t) where {N,TF,SV,SA}

    # extract rate variables
    dλ = dx[SVector{N}(1:N)]
    dvx, dvy, dvz, dωx, dωy, dωz = dx[SVector{6}(N+1:N+6)]
    dδ = dx[end]

    # extract state variables
    λ = x[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = x[SVector{6}(N+1:N+6)]
    δ = x[end]

    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p

    # extract model constants
    bbar = aero.b

    # freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_accelerations(dvx, dvz, dωy)

    # aerodynamic loads
    La, Ma = peters_loads(a, b, ρ, a0, α0, bbar, u, v, ω, vdot, ωdot, λ)

    # loads due to flap deflections
    Lf, Df, Mf = simpleflap_loads(b, u, ρ, clδ, cdδ, cmδ, δ)

    # loads per unit span
    f = SVector(Df, 0, La + Lf)
    m = SVector(0, Ma + Mf, 0)

    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, vdot, ωdot, f..., m...)
end

# --- Performance Overloads --- #

# TODO: state rate, state, and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Peters, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::Peters, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, padd)

    return (rho = padd[1],)
end
