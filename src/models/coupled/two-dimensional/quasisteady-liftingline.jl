
"""
    QuasiSteadyLiftingLine()

Construct a model by coupling a quasi-steady aerodynamic model based on thin airfoil theory 
(see [`QuasiSteady`](@ref)) and a lifting line section model (see [`LiftingLineSection`](@ref)).  
"""
function QuasiSteadyLiftingLine()

    # aerodynamic model
    aero = QuasiSteady()

    # structural model
    stru = LiftingLineSection()

    # submodels
    submodels = (aero, stru)

    # construct coupling
    coupling = quasisteady_liftingline_coupling(aero, stru)

    # return the coupled model
    return CoupledModel(submodels, coupling)
end

# --- Internal Methods for this Coupling --- #

# coupling definition
function quasisteady_liftingline_coupling(aero, stru)

    # coupling function
    g = quasisteady_liftingline_inputs

    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(6)
    ny = Val(6)
    np = Val(8)

    # number of parameters introduced by the coupling (use Val(N) to use inferrable dimensions)
    npc = Val(0)

    # jacobian definitions
    ratejac = Nonlinear() # TODO: define rate jacobian function
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
function quasisteady_liftingline_inputs(dx, x, p, t)
    # rate variables
    dvx, dvy, dvz, dωx, dωy, dωz = dx
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, ρ, c = p
    # freestream velocity components
    u, v, ω = liftinglinesection_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftinglinesection_accelerations(dvx, dvz, dωy)
    # calculate aerodynamic loads
    N, A, M = quasisteady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(f, m)
end