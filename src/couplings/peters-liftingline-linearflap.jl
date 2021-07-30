"""
    couple_models(aero::Peters, stru::LiftingLineSection, flap::LinearFlap,
        ctrl::LiftingLineControl)

Create an aerostructural model using a using the unsteady aerodynamic model
defined by Peters et al, a lifting line aerodynamic model, and a linear steady-state
control surface model.  The existence of this coupling allows [`Peters`](@ref)
and [`LinearFlap`](@ref) to be used with [`LiftingLine`](@ref) and
[`LiftingLineFlaps`](@ref).  This model introduces the freestream air density
``\\rho`` as an additional parameter.
"""
function couple_models(aero::Peters, stru::LiftingLineSection, flap::LinearFlap,
    ctrl::LiftingLineControl)

    return (aero, stru, flap, ctrl)
end

# --- traits --- #

function inplaceness(::Type{Peters}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineControl})

    return OutOfPlace()
end

function mass_matrix_type(::Type{Peters}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineControl})

    return Linear()
end

function state_jacobian_type(::Type{Peters}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineControl})

    return Nonlinear()
end

function number_of_parameters(::Type{Peters}, ::Type{LiftingLineSection},
    ::Type{LinearFlap}, ::Type{LiftingLineControl})

    return 1
end

# --- methods --- #

function get_inputs(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineControl, x, p, t) where {N,TF,SV,SA}
    # extract model constants
    bbar = aero.b
    # extract state variables
    λ = x[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = x[SVector{6}(N+1:N+6)]
    δ = x[end]
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
    # local freestream velocity components
    u = vx
    v = vz
    ω = ωy
    # calculate aerodynamic loads
    L, M = peters_state_loads(a, b, ρ, a0, α0, bbar, u, v, ω, λ)
    # add loads due to flap deflections
    L += ρ*U^2*b*clδ*δ
    D = ρ*U^2*b*cdδ*δ
    M += 2*ρ*U^2*b^2*cmδ*δ
    # forces and moments per unit span
    f = SVector(D, 0, L)
    m = SVector(0, M, 0)
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, 0, 0, f..., m...)
end

function get_input_mass_matrix(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineControl, x, p, t) where {N,TF,SV,SA}

    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p

    # local freestream velocity components
    vdot_dvz = 1
    ωdot_dωy = 1

    # calculate loads
    L_dvx, M_dvx = peters_loads_udot()
    L_dvz, M_dvz = peters_loads_vdot(a, b, ρ)
    L_dωy, M_dωy = peters_loads_ωdot(a, b, ρ)

    # construct submatrices
    Maa = zeros(SMatrix{4,N,TF})
    Mas = @SMatrix [
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -vdot_dvz 0 0 0;
        0 0 0 0 -ωdot_dωy 0]
    Mac = zeros(SMatrix{4,1,TF})

    Msa = zeros(SMatrix{6,N,TF})
    Mss = @SMatrix [
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        -L_dvx 0 -L_dvz 0 -L_dωy 0;
        0 0 0 0 0 0;
        -M_dvx 0 -M_dvz 0 -M_dωy 0;
        0 0 0 0 0 0]
    Msc = zeros(SMatrix{6,1,TF})

    # assemble mass matrix
    return [Maa Mas Mac; Msa Mss Msc]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineControl, x, p, t) where {N,TF,SV,SA}

    # extract model constants
    bbar = aero.b

    # extract state variables
    λ = x[SVector{N}(1:N)]
    vx, vy, vz, ωx, ωy, ωz = x[SVector{6}(N+1:N+6)]
    δ = x[end]

    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p

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

    # add loads due to flap deflections
    L_δ = ρ*U^2*b*clδ
    D_δ = ρ*U^2*b*cdδ
    M_δ += 2*ρ*U^2*b^2*cmδ

    # construct submatrices
    Jaa = zeros(SMatrix{4,N,TF})
    Jas = @SMatrix [u_vx 0 0 0 0 0; 0 0 0 0 ω_ωy 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
    Jac = zeros(SVector{4,1,TF})

    Jsa = vcat(zero(L_λ'), zero(L_λ'), L_λ', zero(M_λ'), M_λ', zero(M_λ'))
    Jss = @SMatrix [0 0 0 0 0 0; 0 0 0 0 0 0; L_vx 0 L_vz 0 L_ωy 0;
        0 0 0 0 0 0; M_vx 0 M_vz 0 M_ωy 0; 0 0 0 0 0 0]
    Jsc = @SVector [D_δ, 0, L_δ, 0, M_δ, 0]

    # assemble jacobian
    return [Jaa Jas Jac; Jsa Jss Jsc]
end

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::Peters{N,TF,SV,SA}, stru::LiftingLineSection,
    flap::LinearFlap, ctrl::LiftingLineControl, dx, x, p, t) where {N,TF,SV,SA}
    # extract state variables
    dλ = dx[SVector{N}(1:N)]
    dvx, dvy, dvz, dωx, dωy, dωz = dx[SVector{6}(N+1:N+6)]
    dδ = dx[end]
    # extract parameters
    a, b, a0, α0, clδ, cdδ, cmδ, ρ = p
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
