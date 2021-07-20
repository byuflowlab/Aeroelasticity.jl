"""
    couple_models(aero::Wagner, stru::LiftingLineSection)

Create an aerostructural model using an unsteady aerodynamic model based
on Wagner's function and a lifting line section model.  The existence of this
coupling allows [`Wagner`](@ref) to be used with [`LiftingLine`](@ref).  This
model introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::Wagner, stru::LiftingLineSection)

# --- traits --- #

inplaceness(::Type{<:Wagner}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Linear()
state_jacobian_type(::Type{<:Wagner}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{<:Wagner}, ::Type{LiftingLineSection}) = 1

# --- methods --- #

function get_inputs(aero::Wagner, stru::LiftingLineSection, s, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # calculate aerodynamic loads
    L, M = wagner_state_loads(a, b, ρ, a0, α0, C1, C2, u, v, ω, λ1, λ2)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return vcat(u, v, ω, f, m)
end

function get_input_mass_matrix(aero::Wagner, stru::LiftingLineSection, s, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
    # calculate loads
    L_dvx, M_dvx = wagner_loads_udot()
    L_dvz, M_dvz = wagner_loads_vdot(a, b, ρ)
    L_dωy, M_dωy = wagner_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0; 0 0; 0 0; 0 0; 0 0]
    Mrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; -L_dvx 0 -L_dvz 0 -L_dωy 0;
        0 0 0 0 0 0; -M_dvx 0 -M_dvz 0 -M_dωy 0; 0 0 0 0 0 0]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Wagner, stru::LiftingLineSection, s, p, t)
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = s
    # extract parameters
    a, b, a0, α0, ρ = p
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
    # compute input jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [u_vx 0 0 0 0 0; 0 0 v_vz 0 0 0; 0 0 0 0 ω_ωy 0]
    Jra = vcat(zero(L_λ'), zero(L_λ'), L_λ', zero(M_λ'), M_λ', zero(M_λ'))
    Jrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_vx 0 L_vz 0 L_ωy 0;
        0 0 0 0 0 0; M_vx 0 M_vz 0 M_ωy 0; 0 0 0 0 0 0]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Wagner, stru::LiftingLineSection,
    ds, s, p, t)
    # extract state rates
    dλ1, dλ2, dvx, dvy, dvz, dωx, dωy, dωz = ds
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    vdot = dvz
    ωdot = dωy
    # calculate aerodynamic loads
    L, M = wagner_rate_loads(a, b, ρ, vdot, ωdot)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(0, 0, 0, f, m)
end
