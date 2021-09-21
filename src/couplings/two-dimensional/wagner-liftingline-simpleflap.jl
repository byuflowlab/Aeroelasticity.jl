"""
    couple_models(aero::Wagner, stru::LiftingLineSection, flap::SimpleFlap,
        ctrl::LiftingLineSectionControl)

Create an aerostructural model using an unsteady aerodynamic model based
on Wagner's function, a lifting line aerodynamic model, and a linear steady-state
control surface model.  The existence of this coupling allows [`Wagner`](@ref)
and [`SimpleFlap`](@ref) to be used with [`LiftingLine`](@ref) and
[`LiftingLineFlaps`](@ref).  This model introduces the freestream air density
``\\rho`` as an additional parameter.
"""
function couple_models(aero::Wagner, stru::LiftingLineSection, flap::SimpleFlap,
    ctrl::LiftingLineSectionControl)

    return (aero, stru, flap, ctrl)
end

# --- Traits --- #

function number_of_additional_parameters(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return 1
end

function coupling_inplaceness(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function coupling_rate_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Linear()
end

function coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_parameter_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function coupling_time_gradient_type(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

# --- Methods --- #

function get_coupling_inputs(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)
    # extract state variables
    dλ1, dλ2, dvx, dvy, dvz, dωx, dωy, dωz, dδ = dx
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, cnδ, caδ, cmδ, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_accelerations(dvx, dvz, dωy)
    # aerodynamic loads
    Na, Aa, Ma = wagner_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # loads due to flap deflections
    Nf, Af, Mf = simpleflap_loads(b, u, ρ, cnδ, caδ, cmδ, δ)
    # loads per unit span
    f = SVector(Aa + Af, 0, Na + Nf)
    m = SVector(0, Ma + Mf, 0)
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, v, ω, f..., m...)
end

# --- Performance Overloads --- #

function get_coupling_rate_jacobian(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)
    # extract state variables
    dλ1, dλ2, dvx, dvy, dvz, dωx, dωy, dωz, dδ = dx
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, cnδ, caδ, cmδ, ρ = p
    # calculate loads
    N_dvx, A_dvx, M_dvx = wagner_loads_udot()
    N_dvz, A_dvz, M_dvz = wagner_loads_vdot(a, b, ρ)
    N_dωy, A_dωy, M_dωy = wagner_loads_ωdot(a, b, ρ)
    # assemble mass matrix
    return @SMatrix [
        0 0      0   0    0   0   0    0   0;
        0 0      0   0    0   0   0    0   0;
        0 0      0   0    0   0   0    0   0;
        0 0    A_dvx 0  A_dvz 0  A_dωy 0   0;
        0 0      0   0    0   0   0    0   0;
        0 0    N_dvx 0  N_dvz 0  N_dωy 0   0;
        0 0      0   0    0   0   0    0   0;
        0 0    M_dvx 0  M_dvz 0  M_dωy 0   0;
        0 0      0   0    0   0   0    0   0
        ]
end

# TODO: State and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, padd)

    return (rho = padd[1],)
end
