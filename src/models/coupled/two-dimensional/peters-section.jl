"""
    peters_section_model(N)

Construct a model by coupling an unsteady aerodynamic model based on Peters' finite state 
theory (see [`Peters`](@ref)) and a two-degree of freedom typical section model (see 
[`Section`](@ref)).  This model introduces the freestream velocity ``U_\\infty``, 
air density ``\\rho_\\infty`` and air speed of sound ``c`` as additional parameters.
"""
function peters_section_model(N)

    # aerodynamic model
    aero = peters_model(N)

    # structural model
    stru = typical_section_model()

    # submodels
    submodels = (aero, stru)

    # construct coupling
    coupling = peters_section_coupling(aero, stru)

    # return the coupled model
    return CoupledModel(submodels, coupling)
end

# --- Internal Methods for this Coupling --- #

# coupling definition
function peters_section_coupling(aero, stru)

    # state variable indices
    iλ, iq = state_indices((aero, stru))

    # model constants
    bbar = aero.constants.bbar

    # number of aerodynamic state variables
    N = length(bbar)

    # coupling function
    g = (dx, x, p, t) -> peters_section_inputs(dx, x, p, t; iλ, iq, bbar)

    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(N+4)
    ny = Val(6)
    np = Val(14)

    # number of parameters introduced by the coupling (use Val(N) to use inferrable dimensions)
    npc = Val(3)

    # jacobian definitions
    ratejac = Nonlinear() # TODO: define rate jacobian function
    statejac = Nonlinear() # TODO: define state jacobian function
    paramjac = Nonlinear() # TODO: define parameter jacobian function
    tgrad = Zeros()
    
    # convenience function for setting coupling parameters
    setparam = peters_section_setparam

    # convenience function for separating coupling parameters
    sepparam = peters_section_sepparam
    
    # return resulting coupling
    return Coupling{false}(g, nx, ny, np, npc;
        ratejac = ratejac,
        statejac = statejac,
        paramjac = paramjac,
        tgrad = tgrad,
        setparam = setparam,
        sepparam = sepparam)
end

# coupling function
function peters_section_inputs(dx, x, p, t; iλ, iq, bbar)
    # extract rate variables
    dλ = view(dx, iλ)
    dh, dθ, dhdot, dθdot = view(dx, iq)
    # extract state variables
    λ = view(x, iλ)
    h, θ, hdot, θdot = view(x, iq)
    # extract parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ, c = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    N, A, M = peters_loads(a, b, ρ, c, a0, α0, cd0, cm0, bbar, u, v, ω, vdot, ωdot, λ)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return SVector(u, ω, vdot, ωdot, L, M)
end

# convenience function for defining the coupling function parameters
function peters_section_setparam(p; U, rho, c)
    p[1] = U
    p[2] = rho
    p[3] = c
    return p
end

# convenience function for separating the coupling function parameters
peters_section_sepparam(p) = (U = p[1], rho = p[2], c = p[3])
