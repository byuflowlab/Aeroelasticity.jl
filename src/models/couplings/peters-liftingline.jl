
"""
    Coupling(models::Tuple{Peters, LiftingLineSection})

Construct a model by coupling an unsteady aerodynamic model based on Peters' finite state 
theory (see [`Peters`](@ref)) and a lifting line section model (see 
[`LiftingLineSection`](@ref)).
"""
function Coupling(models::Tuple{Peters{N}, LiftingLineSection}, submodels=Submodel.(models)) where N

    peters = models[1]

    # state variable indices
    iλ = 1:N
    iq = N+1:N+6

    # coupling function
    g = (dx, x, p, t) -> peters_liftingline_inputs(dx, x, p, t; iλ, iq, bbar = peters.b)

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
    u, v, ω = liftingline_section_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_section_accelerations(dvx, dvz, dωy)
    # calculate aerodynamic loads
    N, A, M = peters_loads(a, b, ρ, c, a0, α0, cd0, cm0, bbar, u, v, ω, vdot, ωdot, λ)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return SVector(u, ω, vdot, ωdot, f..., m...)
end