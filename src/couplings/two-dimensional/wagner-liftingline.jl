"""
    couple_models(aero::Wagner, stru::LiftingLineSection)

Create an aerostructural model using an unsteady aerodynamic model based
on Wagner's function and a lifting line section model.  The existence of this
coupling allows [`Wagner`](@ref) to be used with [`LiftingLine`](@ref).  This
model introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::Wagner, stru::LiftingLineSection) = (aero, stru)

# --- Traits --- #

number_of_additional_parameters(::Type{<:Wagner}, ::Type{LiftingLineSection}) = 1

coupling_inplaceness(::Type{<:Wagner}, ::Type{LiftingLineSection}) = OutOfPlace()

coupling_rate_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Linear()
coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Nonlinear()
coupling_time_gradient_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Zeros()

# --- Methods --- #

function get_coupling_inputs(aero::Wagner, stru::LiftingLineSection, dx, x, p, t)
    # extract rate variables
    dλ1, dλ2, dvx, dvy, dvz, dωx, dωy, dωz = dx
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = x
    # extract parameters
    a, b, a0, α0, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_accelerations(dvx, dvz, dωy)
    # aerodynamic loads
    N, A, M = wagner_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return vcat(u, v, ω, f, m)
end

# --- Performance Overloads --- #

# TODO: State rate, state, and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::Wagner, stru::LiftingLineSection; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::Wagner, stru::LiftingLineSection, padd)

    return (rho = padd[1],)
end
