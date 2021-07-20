"""
    couple_models(aero::Peters, stru::LiftingLineSection)

Create an aerostructural model using a using the unsteady aerodynamic model
defined by Peters et al. and a lifting line section model.  The existence of this
coupling allows [`Peters`](@ref) to be used with [`LiftingLine`](@ref).  This
model introduces the freestream air density ``\\rho`` as an additional parameter.
"""
couple_models(aero::Peters, stru::LiftingLineSection)

# --- traits --- #

inplaceness(::Type{<:Peters}, ::Type{LiftingLineSection}) = OutOfPlace()
mass_matrix_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Linear()
state_jacobian_type(::Type{<:Peters}, ::Type{LiftingLineSection}) = Nonlinear()
number_of_parameters(::Type{<:Peters}, ::Type{LiftingLineSection}) = 1

# --- methods --- #

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = s[SVector{6}(N+1:N+6)]
    # extract parameters
    a, b, a0, α0, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = peters_state_loads(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return vcat(u, ω, 0, 0, f, m)
end

function get_input_mass_matrix(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    vdot_dvz = 1
    ωdot_dωy = 1
    # calculate loads
    L_dvx, M_dvx = peters_loads_udot()
    L_dvz, M_dvz = peters_loads_vdot(a, b, ρ)
    L_dωy, M_dωy = peters_loads_ωdot(a, b, ρ)
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0;
        0 0 -vdot_dvz 0 0 0; 0 0 0 0 -ωdot_dωy 0]
    Mra = zeros(SMatrix{6,N,TF})
    Mrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; -L_dvx 0 -L_dvz 0 -L_dωy 0;
        0 0 0 0 0 0; -M_dvx 0 -M_dvz 0 -M_dωy 0; 0 0 0 0 0 0]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = s[SVector{6}(N+1:N+6)]
    # extract parameters
    a, b, a0, α0, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    u_vx = 1
    ω_ωy = 1
    # compute loads
    out = peters_loads_λ(a, b, ρ, a0, bbar, u)
    L_λ, M_λ = out[1,:], out[2,:]
    L_vx, M_vx = peters_loads_u(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    L_vz, M_vz = peters_loads_v(a, b, ρ, a0, u)
    L_ωy, M_ωy = peters_loads_ω(a, b, ρ, a0, u)
    # construct submatrices
    Jda = zeros(SMatrix{4,N,TF})
    Jds = @SMatrix [u_vx 0 0 0 0 0; 0 0 0 0 ω_ωy 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
    Jra = vcat(zero(L_λ'), zero(L_λ'), L_λ', zero(M_λ'), M_λ', zero(M_λ'))
    Jrs = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_vx 0 L_vz 0 L_ωy 0;
        0 0 0 0 0 0; M_vx 0 M_vz 0 M_ωy 0; 0 0 0 0 0 0]
    # assemble jacobian
    return [Jda Jds; Jra Jrs]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Peters{N,TF,SV,SA},
    stru::LiftingLineSection, ds, s, p, t) where {N,TF,SV,SA}
    # extract state rates
    dλ = ds[SVector{N}(1:N)]
    dvx, dvy, dvz, dωx, dωy, dωz = ds[SVector{6}(N+1:N+6)]
    # extract parameters
    a, b, a0, α0, ρ = p
    # local freestream velocity components
    vdot = dvz
    ωdot = dωy
    # calculate aerodynamic loads
    L, M = peters_rate_loads(a, b, ρ, vdot, ωdot)
    # forces and moments per unit span
    f = SVector(0, 0, L)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(0, 0, vdot, ωdot, f, m)
end
