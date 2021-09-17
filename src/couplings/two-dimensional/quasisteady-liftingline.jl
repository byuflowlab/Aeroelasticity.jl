"""
    couple_models(aero::QuasiSteady, stru::LiftingLineSection)

Create an aerostructural model using a quasi-steady aerodynamics model and a
lifting line section model.  The existence of this coupling allows
[`QuasiSteady`](@ref) to be used with [`LiftingLine`](@ref). This model
introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::QuasiSteady, stru::LiftingLineSection) = (aero, stru)

# --- Traits --- #

# steady
number_of_additional_parameters(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = 1
coupling_inplaceness(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Nonlinear()
copuling_time_gradient_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Zeros()

# quasisteady without acceleration terms
number_of_additional_parameters(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = 1
coupling_inplaceness(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Nonlinear()
copuling_time_gradient_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Zeros()

# quasisteady with acceleration terms
number_of_additional_parameters(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = 1
coupling_inplaceness(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = OutOfPlace()
coupling_rate_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Linear()
coupling_state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Nonlinear()
coupling_parameter_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Nonlinear()
copuling_time_gradient_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Zeros()

# --- Methods --- #

# steady
function get_coupling_inputs(aero::QuasiSteady{0}, stru::LiftingLineSection, dx, x, p, t)
    # extract rate variables
    dvx, dvy, dvz, dωx, dωy, dωz = dx
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = x
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    # calculate loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # loads per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return SVector(f..., m...)
end

# quasisteady without accelaration terms
function get_coupling_inputs(aero::QuasiSteady{1}, stru::LiftingLineSection, dx, x, p, t)
    # extract rate variables
    vx, vy, vz, ωx, ωy, ωz = dx
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = x
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # loads per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return SVector(f..., m...)
end

# quasisteady with acceleration terms
function get_coupling_inputs(aero::QuasiSteady{2}, stru::LiftingLineSection, dx, x, p, t)
    # extract rate variables
    dvx, dvy, dvz, dωx, dωy, dωz = dx
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = x
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u, v, ω = liftingline_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_accelerations(dvx, dvz, dωy)
    # calculate aerodynamic loads
    L, M = quasisteady2_loads(a, b, ρ, a0, α0, u, v, ω, vdot, ωdot)
    # loads per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return SVector(f..., m...)
end

# --- Performance Overloads --- #

# TODO: State rate, state, and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, aero::QuasiSteady, stru::LiftingLineSection; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::QuasiSteady, stru::LiftingLineSection, padd)

    return (rho = padd[1],)
end
