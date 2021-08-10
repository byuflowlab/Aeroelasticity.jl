"""
    couple_models(aero::QuasiSteady, stru::LiftingLineSection, flap::LinearFlap,
        ctrl::LiftingLineSectionControl)

Create an aerostructural model using a quasi-steady aerodynamics model, a
lifting line aerodynamic model, and a linear, steady-state control surface model.
The existence of this coupling allows the [`QuasiSteady`](@ref) and [`LinearFlap`](@ref)
to be used with [`LiftingLine`](@ref) and [`LiftingLineFlaps`](@ref). This model
introduces the freestream air density ``\\rho`` as an additional parameter.
"""
function couple_models(aero::QuasiSteady, stru::LiftingLineSection, flap::LinearFlap,
    ctrl::LiftingLineSectionControl)

    return (aero, stru, flap, ctrl)
end

# --- traits --- #

function inplaceness(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function mass_matrix_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

function state_jacobian_type(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function number_of_parameters(::Type{QuasiSteady{0}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return 1
end

function inplaceness(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function mass_matrix_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return Zeros()
end

function state_jacobian_type(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function number_of_parameters(::Type{QuasiSteady{1}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return 1
end

function inplaceness(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function mass_matrix_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return Linear()
end

function state_jacobian_type(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

function number_of_parameters(::Type{QuasiSteady{2}}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineSectionControl})
    return 1
end

# --- methods --- #

function get_inputs(aero::QuasiSteady{0}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    # calculate loads
    L, M = quasisteady0_loads(a, b, ρ, a0, α0, u, v)
    # add loads due to flap deflections
    L += ρ*u^2*b*clδ*δ
    D = ρ*u^2*b*cdδ*δ
    M += 2*ρ*u^2*b^2*cmδ*δ
    # forces and moments per unit span
    f = SVector(D, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_inputs(aero::QuasiSteady{1}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = quasisteady1_loads(a, b, ρ, a0, α0, u, v, ω)
    # add loads due to flap deflections
    L += ρ*u^2*b*clδ*δ
    D = ρ*u^2*b*cdδ*δ
    M += 2*ρ*u^2*b^2*cmδ*δ
    # forces and moments per unit span
    f = SVector(D, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_inputs(aero::QuasiSteady{2}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = quasisteady2_state_loads(a, b, ρ, a0, α0, u, v, ω)
    # add loads due to flap deflections
    L += ρ*u^2*b*clδ*δ
    D = ρ*u^2*b*cdδ*δ
    M += 2*ρ*u^2*b^2*cmδ*δ
    # forces and moments per unit span
    f = SVector(D, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end

function get_input_mass_matrix(aero::QuasiSteady{2}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # calculate loads
    L_vx, M_vx = quasisteady2_udot()
    L_vz, M_vz = quasisteady2_vdot(a, b, ρ)
    L_ωy, M_ωy = quasisteady2_ωdot(a, b, ρ)
    # return input mass matrix
    return @SMatrix [0 0 0 0 0 0 0; 0 0 0 0 0 0 0; -L_vx 0 -L_vz 0 -L_ωy 0 0; 0 0 0 0 0 0 0;
        -M_vx 0 -M_vz 0 -M_ωy 0 0; 0 0 0 0 0 0 0]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::QuasiSteady{0}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    # calculate loads
    L_u, M_u = quasisteady0_u(a, b, ρ, a0, v)
    L_v, M_v = quasisteady0_v(a, b, ρ, a0, u)
    # add loads due to flap deflections
    L_u += 2*ρ*u*b*clδ*δ
    D_u = 2*ρ*u*b*cdδ*δ
    M_u += 4*ρ*u*b^2*cmδ*δ
    L_δ = ρ*u^2*b*clδ
    D_δ = ρ*u^2*b*cdδ
    M_δ = 2*ρ*u^2*b^2*cmδ
    # return inputs
    return @SMatrix [D_u 0 0 0 0 0 D_δ; 0 0 0 0 0 0 0; L_u 0 L_v 0 0 0 L_δ; 0 0 0 0 0 0 0;
        M_u 0 M_v 0 0 0 M_δ; 0 0 0 0 0 0 0]
end

function get_input_state_jacobian(aero::QuasiSteady{1}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate loads
    L_u, M_u = quasisteady1_u(a, b, ρ, a0, α0, u, v, ω)
    L_v, M_v = quasisteady1_v(a, b, ρ, a0, α0, u)
    L_ω, M_ω = quasisteady1_ω(a, b, ρ, a0, u)
    # add loads due to flap deflections
    L_u += 2*ρ*u*b*clδ*δ
    D_u = 2*ρ*u*b*cdδ*δ
    M_u += 4*ρ*u*b^2*cmδ*δ
    L_δ = ρ*u^2*b*clδ
    D_δ = ρ*u^2*b*cdδ
    M_δ = 2*ρ*u^2*b^2*cmδ
    # return inputs
    return @SMatrix [D_u 0 0 0 0 0 D_δ; 0 0 0 0 0 0 0; L_u 0 L_v 0 L_ω 0 L_δ;
        0 0 0 0 0 0 0; M_u 0 M_v 0 M_ω 0 M_δ; 0 0 0 0 0 0 0]
end

function get_input_state_jacobian(aero::QuasiSteady{2}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate loads
    L_u, M_u = quasisteady2_u(a, b, ρ, a0, α0, u, v, ω)
    L_v, M_v = quasisteady2_v(a, b, ρ, a0, α0, u)
    L_ω, M_ω = quasisteady2_ω(a, b, ρ, a0, u)
    # add loads due to flap deflections
    L_u += 2*ρ*u*b*clδ*δ
    D_u = 2*ρ*u*b*cdδ*δ
    M_u += 4*ρ*u*b^2*cmδ*δ
    L_δ = ρ*u^2*b*clδ
    D_δ = ρ*u^2*b*cdδ
    M_δ = 2*ρ*u^2*b^2*cmδ
    # return inputs
    return @SMatrix [D_u 0 0 0 0 0 D_δ; 0 0 0 0 0 0 0; L_u 0 L_v 0 L_ω 0 L_δ;
        0 0 0 0 0 0 0; M_u 0 M_v 0 M_ω 0 M_δ; 0 0 0 0 0 0 0]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::QuasiSteady{0}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)

    return @SVector zeros(6)
end

function get_inputs_from_state_rates(aero::QuasiSteady{1}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)

    return @SVector zeros(6)
end

function get_inputs_from_state_rates(aero::QuasiSteady{2}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)
    # extract state rates
    dvx, dvy, dvz, dωx, dωy, dωz, dδ = dx
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
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
