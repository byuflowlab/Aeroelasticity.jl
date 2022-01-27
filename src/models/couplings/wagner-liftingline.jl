
"""
    wagner_liftingline_model(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

Construct a model by coupling an unsteady aerodynamic model based on Wagner's function 
(see [`wagner_model`](@ref)) and a lifting line section model (see 
[`liftingline_section_model`](@ref)).
"""
function wagner_liftingline_model(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

    # aerodynamic model
    aero = wagner_model(; C1, C2, eps1, eps2)

    # structural model
    stru = liftingline_section_model()

    # submodels
    submodels = (aero, stru)

    # construct coupling
    coupling = wagner_liftingline_coupling(aero, stru)

    # return the coupled model
    return CoupledModel(submodels, coupling)
end

# --- Internal Methods for this Coupling --- #

# coupling definition
function wagner_liftingline_coupling(aero, stru)

    # model constants
    C1 = aero.constants.C1
    C2 = aero.constants.C2

    # coupling function
    g = (dx, x, p, t) -> wagner_liftingline_inputs(dx, x, p, t; C1, C2)
    
    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(8)
    ny = Val(9)
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
function wagner_liftingline_inputs(dx, x, p, t; C1, C2)
    # rate variables
    dλ1, dλ2, dvx, dvy, dvz, dωx, dωy, dωz = dx
    # extract state variables
    λ1, λ2, vx, vy, vz, ωx, ωy, ωz = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, ρ, c = p
    # freestream velocity components
    u, v, ω = liftingline_section_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_section_accelerations(dvx, dvz, dωy)
    # calculate aerodynamic loads
    N, A, M = wagner_loads(a, b, ρ, c, a0, α0, cd0, cm0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return vcat(u, v, ω, f, m)
end