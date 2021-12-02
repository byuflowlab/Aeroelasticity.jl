
"""
    PetersLiftingLine()

Construct a model by coupling an unsteady aerodynamic model based on Peters' finite state 
theory (see [`Peters`](@ref)) and a lifting line section model (see 
[`LiftingLineSection`](@ref)).
"""
function PetersLiftingLine(N)

    # aerodynamic model
    aero = Peters(N)

    # structural model
    stru = LiftingLineSection()

    # submodels
    submodels = (aero, stru)

    # construct coupling
    coupling = peters_liftingline_coupling(aero, stru)

    # return the coupled model
    return CoupledModel(submodels, coupling)
end

# --- Internal Methods for this Coupling --- #

# coupling definition
function peters_liftingline_coupling(aero, stru)

    # state variable indices
    iλ, iq = state_indices((aero, stru))

    # model constants
    bbar = aero.constants.bbar

    # number of aerodynamic state variables
    N = length(bbar)

    # coupling function
    g = (dx, x, p, t) -> peters_liftingline_inputs(dx, x, p, t; iλ, iq, bbar)

    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(N+6)
    ny = Val(10)
    np = Val(8)

    # number of parameters introduced by the coupling (use Val(N) to use inferrable dimensions)
    npc = Val(0)

    # jacobian definitions
    ratejac = Linear() # TODO: define rate jacobian function
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
function peters_liftingline_inputs(dx, x, p, t; iλ, iq, bbar)
    # extract rate variables
    dλ = view(dx, iλ)
    dvx, dvy, dvz, dωx, dωy, dωz = view(dx, iq)
    # extract state variables
    λ = view(x, iλ)
    vx, vy, vz, ωx, ωy, ωz = view(x, iq)
    # extract parameters
    a, b, a0, α0, cd0, cm0, ρ, c = p
    # freestream velocity components
    u, v, ω = liftinglinesection_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftinglinesection_accelerations(dvx, dvz, dωy)
    # calculate aerodynamic loads
    N, A, M = peters_loads(a, b, ρ, c, a0, α0, cd0, cm0, bbar, u, v, ω, vdot, ωdot, λ)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return SVector(u, ω, vdot, ωdot, f..., m...)
end