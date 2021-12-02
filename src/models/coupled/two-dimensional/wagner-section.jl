"""
    WagnerSection(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)

Construct a model by coupling an unsteady aerodynamic model based on Wagner's function (see 
[`Wagner`](@ref)) and a two-degree of freedom typical section model (see [`Section`](@ref)).  
This coupling introduces the freestream velocity ``U_\\infty``, air density 
``\\rho_\\infty`` and air speed of sound ``c`` as additional parameters.
"""
function WagnerSection(; C1=0.165, C2=0.335, eps1 = 0.0455, eps2 = 0.3)
    
    # aerodynamic model
    aero = Wagner(; C1, C2, eps1, eps2)

    # structural model
    stru = Section()

    # submodels
    submodels = (aero, stru)

    # construct coupling
    coupling = wagner_section_coupling(aero, stru)

    # return the coupled model
    return CoupledModel(submodels, coupling)
end

# --- Internal Methods for this Coupling --- #

# coupling definition
function wagner_section_coupling(aero, stru)

    # model constants
    C1 = aero.constants.C1
    C2 = aero.constants.C2

    # coupling function
    g = (dx, x, p, t) -> wagner_section_inputs(dx, x, p, t; C1, C2)
    
    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(6)
    ny = Val(5)
    np = Val(14)

    # number of parameters introduced by the coupling (use Val(N) to use inferrable dimensions)
    npc = Val(3)

    # jacobian definitions
    ratejac = Nonlinear() # TODO: define rate jacobian function
    statejac = Nonlinear() # TODO: define state jacobian function
    paramjac = Nonlinear() # TODO: define parameter jacobian function
    tgrad = Zeros()

    # convenience function for setting coupling parameters
    setparam = wagner_section_setparam

    # convenience function for separating coupling parameters
    sepparam = wagner_section_sepparam

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
function wagner_section_inputs(dx, x, p, t; C1, C2)
    # extract rate variables
    dλ1, dλ2, dh, dθ, dhdot, dθdot = dx
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = x
    # extract parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ, c = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    N, A, M = wagner_loads(a, b, ρ, c, a0, α0, cd0, cm0, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # lift is approximately normal force
    L = N
    # return inputs
    return SVector(u, v, ω, L, M)
end

# convenience function for defining the coupling function parameters
function wagner_section_setparam(p; U, rho, c)
    p[1] = U
    p[2] = rho
    p[3] = c
    return p
end

# convenience function for separating the coupling function parameters
wagner_section_sepparam(p) = (U = p[1], rho = p[2], c = p[3])
