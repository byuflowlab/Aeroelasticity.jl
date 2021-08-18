"""
    couple_models(aero::QuasiSteady, stru::LiftingLineSection)

Create an aerostructural model using a quasi-steady aerodynamics model and a
lifting line section model.  The existence of this coupling allows
[`QuasiSteady`](@ref) to be used with [`LiftingLine`](@ref). This model
introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::QuasiSteady, stru::LiftingLineSection) = (aero, stru)

# --- traits --- #

number_of_additional_parameters(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = 1
coupling_inplaceness(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = OutOfPlace()
coupling_mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection}) = Nonlinear()

number_of_additional_parameters(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = 1
coupling_inplaceness(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = OutOfPlace()
coupling_mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Zeros()
coupling_state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection}) = Nonlinear()

number_of_additional_parameters(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = 1
coupling_inplaceness(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = OutOfPlace()
coupling_mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Linear()
coupling_state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection}) = Nonlinear()

# --- methods --- #

function get_coupling_inputs(aero::QuasiSteady{0}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    # calculate loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_coupling_inputs(aero::QuasiSteady{1}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_coupling_inputs(aero::QuasiSteady{2}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, ω)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_coupling_mass_matrix(aero::QuasiSteady{2}, stru::LiftingLineSection,
    s, p, t)
    # extract parameters
    a, b, a0, α0, ρ = p
    # calculate loads
    L_vx, M_vx = quasisteady2_udot()
    L_vz, M_vz = quasisteady2_vdot(a, b, ρ)
    L_ωy, M_ωy = quasisteady2_ωdot(a, b, ρ)
    # return input mass matrix
    return @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; -L_vx 0 -L_vz 0 -L_ωy 0; 0 0 0 0 0 0;
        -M_vx 0 -M_vz 0 -M_ωy 0; 0 0 0 0 0 0]
end

# --- performance overloads --- #

function get_coupling_state_jacobian(aero::QuasiSteady{0}, stru::LiftingLineSection, s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    # calculate loads
    L_u, M_u = quasisteady0_u(a, b, ρ, a0, v)
    L_v, M_v = quasisteady0_v(a, b, ρ, a0, u)
    # return inputs
    return @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_u 0 L_v 0 0 0; 0 0 0 0 0 0;
        M_u 0 M_v 0 0 0; 0 0 0 0 0 0]
end

function get_coupling_state_jacobian(aero::QuasiSteady{1}, stru::LiftingLineSection,
    s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate loads
    L_u, M_u = quasisteady1_u(a, b, ρ, a0, α0, u, v, ω)
    L_v, M_v = quasisteady1_v(a, b, ρ, a0, α0, u)
    L_ω, M_ω = quasisteady1_ω(a, b, ρ, a0, u)
    # return inputs
    return @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_u 0 L_v 0 L_ω 0;
        0 0 0 0 0 0; M_u 0 M_v 0 M_ω 0; 0 0 0 0 0 0]
end

function get_coupling_state_jacobian(aero::QuasiSteady{2}, stru::LiftingLineSection,
    s, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate loads
    L_u, M_u = quasisteady2_u(a, b, ρ, a0, α0, u, v, ω)
    L_v, M_v = quasisteady2_v(a, b, ρ, a0, α0, u)
    L_ω, M_ω = quasisteady2_ω(a, b, ρ, a0, u)
    # return inputs
    return @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_u 0 L_v 0 L_ω 0; 0 0 0 0 0 0;
        M_u 0 M_v 0 M_ω 0; 0 0 0 0 0 0]
end

# --- unit testing methods --- #

function get_coupling_inputs_using_state_rates(aero::QuasiSteady{0}, stru::LiftingLineSection,
    ds, s, p, t)

    return @SVector zeros(6)
end

function get_coupling_inputs_using_state_rates(aero::QuasiSteady{1}, stru::LiftingLineSection,
    ds, s, p, t)

    return @SVector zeros(6)
end

function get_coupling_inputs_using_state_rates(aero::QuasiSteady{2}, stru::LiftingLineSection,
    ds, s, p, t)
    # extract state rates
    dvx, dvy, dvz, dωx, dωy, dωz = ds
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    du = dvx
    dv = dvz
    dω = dωy
    # calculate aerodynamic loads
    L, M = quasisteady2_rate_loads(a, b, ρ, dv, dω)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

# --- convenience methods --- #

function set_additional_parameters!(padd, aero::QuasiSteady, stru::LiftingLineSection; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::QuasiSteady, stru::LiftingLineSection, padd)

    return (rho = padd[1],)
end
