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

# --- traits --- #

function number_of_additional_parameters(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return 1
end

function coupling_inplaceness(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return OutOfPlace()
end

function coupling_mass_matrix_type(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Linear()
end

function coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection},
    ::Type{SimpleFlap}, ::Type{LiftingLineSectionControl})

    return Nonlinear()
end

# --- methods --- #

function get_inputs(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = wagner_state_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # add loads due to flap deflections
    L += ρ*u^2*b*clδ*δ
    D = ρ*u^2*b*cdδ*δ
    M += 2*ρ*u^2*b^2*cmδ*δ
    # forces and moments per unit span
    f = SVector(D, 0, L)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, v, ω, f..., m...)
end

function get_coupling_mass_matrix(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # calculate loads
    L_dvx, M_dvx = wagner_loads_udot()
    L_dvz, M_dvz = wagner_loads_vdot(a, b, ρ)
    L_dωy, M_dωy = wagner_loads_ωdot(a, b, ρ)
    # assemble mass matrix
    return @SMatrix [
        0 0      0   0    0   0   0    0   0;
        0 0      0   0    0   0   0    0   0;
        0 0      0   0    0   0   0    0   0;
        0 0      0   0    0   0   0    0   0;
        0 0      0   0    0   0   0    0   0;
        0 0   -L_dvx 0 -L_dvz 0 -L_dωy 0   0;
        0 0      0   0    0   0   0    0   0;
        0 0   -M_dvx 0 -M_dvz 0 -M_dωy 0   0;
        0 0      0   0    0   0   0    0   0
        ]
end

# --- performance overloads --- #

function get_coupling_state_jacobian(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, x, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz, δ = x
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    u_vx = 1
    v_vz = 1
    ω_ωy = 1
    # model constants
    C1 = aero.C1
    C2 = aero.C2
    # calculate loads
    out = wagner_loads_λ(a, b, ρ, a0, u)
    L_λ, M_λ = out[1,:], out[2,:]
    L_vx, M_vx = wagner_loads_u(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    L_vz, M_vz = wagner_loads_v(a, b, ρ, a0, C1, C2, u)
    L_ωy, M_ωy = wagner_loads_ω(a, b, ρ, a0, C1, C2, u)
    # add loads due to flap deflections
    L_vx += 2*ρ*u*b*clδ*δ
    D_vx = 2*ρ*u*b*cdδ*δ
    M_vx += 4*ρ*u*b^2*cmδ*δ
    L_δ = ρ*u^2*b*clδ
    D_δ = ρ*u^2*b*cdδ
    M_δ = 2*ρ*u^2*b^2*cmδ
    # assemble jacobian
    return @SMatrix [
         0      0     u_vx 0   0   0   0    0  0;
         0      0      0   0 v_vz  0   0    0  0;
         0      0      0   0   0   0  ω_ωy  0  0;

         0      0     D_vx 0   0   0   0    0 D_δ;
         0      0      0   0   0   0   0    0  0;
        L_λ[1] L_λ[2] L_vx 0 L_vz  0  L_ωy  0 L_δ;
         0      0      0   0   0   0   0    0  0;
        M_λ[1] M_λ[2] M_vx 0 M_vz  0  M_ωy  0 M_δ;
         0      0      0   0   0   0   0    0  0
        ]
end

# --- unit testing methods --- #

function get_inputs_using_state_rates(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, dx, x, p, t)
    # extract state rates
    dλ1, dλ2, dvx, dvy, dvz, dωx, dωy, dωz, dδ = dx
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    vdot = dvz
    ωdot = dωy
    # calculate aerodynamic loads
    L, M = wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return SVector(0, 0, 0, f..., m...)
end

# --- convenience methods --- #

function set_additional_parameters!(padd, aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl; rho)

    padd[1] = rho

    return padd
end

function separate_additional_parameters(aero::Wagner, stru::LiftingLineSection,
    flap::SimpleFlap, ctrl::LiftingLineSectionControl, padd)

    return (rho = padd[1],)
end
