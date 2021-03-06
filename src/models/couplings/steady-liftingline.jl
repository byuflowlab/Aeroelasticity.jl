"""
    steady_liftingline_model()

Construct a model by coupling a steady aerodynamic model based on thin airfoil theory (see 
[`Steady`](@ref)) and a lifting line section model (see 
[`liftingline_section_model`](@ref)).  
"""
function steady_liftingline_model()

    # aerodynamic model
    aero = steady_model()

    # structural model
    stru = liftingline_section_model()

    # submodels
    submodels = (aero, stru)

    # construct coupling
    coupling = steady_liftingline_coupling(aero, stru)

    # return the coupled model
    return CoupledModel(submodels, coupling)
end

# --- Internal Methods --- #

function steady_liftingline_coupling(aero, stru)

    # coupling function
    g = steady_liftingline_inputs

    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(6)
    ny = Val(6)
    np = Val(8)

    # number of parameters introduced by the coupling (use Val(N) to use inferrable dimensions)
    npc = Val(0)

    # jacobian definitions
    ratejac = Zeros()
    statejac = Nonlinear() # TODO: define state jacobian function
    paramjac = Nonlinear() # TODO: define parameter jacobian function
    tgrad = Zeros()

    # return resulting coupling
    return Coupling{false}(g, nx, ny, np, npc;
        ratejac = ratejac,
        statejac = statejac,
        paramjac = paramjac,
        tgrad = tgrad)
end

# coupling function
function steady_liftingline_inputs(dx, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, ρ, c = p
    # freestream velocity components
    u, v, ω = liftingline_section_velocities(vx, vz, ωy)
    # calculate aerodynamic loads
    N, A, M = steady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end