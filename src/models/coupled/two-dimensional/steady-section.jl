"""
    steady_section_model()

Construct a model by coupling a steady aerodynamic model based on thin airfoil theory (see 
[`Steady`](@ref)) and a two-degree of freedom typical section model (see [`Section()`]).  
This model introduces the freestream velocity ``U_\\infty``, air density ``\\rho_\\infty``, 
and air speed of sound ``c`` as additional parameters.
"""
function steady_section_model()

    # aerodynamic model
    aero = steady_model()

    # structural model
    stru = typical_section_model()

    # submodels
    submodels = (aero, stru)

    # construct coupling
    coupling = steady_section_coupling(aero, stru)

    # return the coupled model
    return CoupledModel(submodels, coupling)
end

# --- Internal Methods for this Coupling --- #

# coupling definition
function steady_section_coupling(aero, stru)

    # coupling function
    g = steady_section_inputs

    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(4)
    ny = Val(2)
    np = Val(14)

    # number of parameters introduced by the coupling (use Val(N) to use inferrable dimensions)
    npc = Val(3)

    # jacobian definitions
    ratejac = Zeros()
    statejac = Nonlinear() # TODO: define state jacobian function
    paramjac = Nonlinear() # TODO: define parameter jacobian function
    tgrad = Zeros()

    # convenience function for setting coupling parameters
    setparam = steady_section_setparam

    # convenience function for separating coupling parameters
    sepparam = steady_section_sepparam

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
function steady_section_inputs(dx, x, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ, c = p
    # local freestream velocity components
    u, v = section_steady_velocities(U, θ)
    # calculate aerodynamic loads
    N, A, M = steady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return SVector(L, M)
end

# convenience function for defining the coupling function parameters
function steady_section_setparam(p; U, rho, c)
    p[1] = U
    p[2] = rho
    p[3] = c
    return p
end

# convenience function for separating the coupling function parameters
steady_section_sepparam(p) = (U = p[1], rho = p[2], c = p[3])

# local freestream velocities
function section_steady_velocities(U, θ)
    u = U
    v = U*θ
    return SVector(u, v)
end