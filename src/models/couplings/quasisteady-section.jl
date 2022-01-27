# --- Coupling Model Creation --- #

"""
    Coupling(::QuasiSteady, ::Section)

Coupling model for coupling a quasi-steady aerodynamic model based on thin airfoil theory 
(see [`QuasiSteady`](@ref)) and a two-degree of freedom typical section model 
(see [`Section()`]).  This model introduces the freestream velocity ``U``, air density 
``\\rho``, and air speed of sound ``c`` as additional parameters.
"""
function Coupling(::QuasiSteady, ::Section)

    # coupling function
    g = quasisteady_section_inputs

    # number of states, inputs, and parameters (use Val(N) to use inferrable dimensions)
    nx = Val(4)
    ny = Val(2)
    np = Val(14)

    # number of parameters introduced by the coupling (use Val(N) to use inferrable dimensions)
    npc = Val(3)

    # jacobian definitions
    ratejac = Nonlinear() # TODO: define rate jacobian function
    statejac = Nonlinear() # TODO: define state jacobian function
    paramjac = Nonlinear() # TODO: define parameter jacobian function
    tgrad = Zeros()

    # convenience function for setting coupling parameters
    setparam = quasisteady_section_set_parameters!

    # convenience function for separating coupling parameters
    sepparam = quasisteady_section_separate_parameters

    # return resulting coupling
    return Coupling{false}(g, nx, ny, np, npc;
        ratejac = ratejac,
        statejac = statejac,
        paramjac = paramjac,
        tgrad = tgrad,
        setparam = setparam,
        sepparam = sepparam)
end

# --- Internal Methods --- #

# coupling function
function quasisteady_section_inputs(dx, x, p, t)
    # extract state variables
    dh, dθ, dhdot, dθdot = dx
    # extract state variables
    h, θ, hdot, θdot = x
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, cd0, cm0, kh, kθ, m, Sθ, Iθ, U, ρ, c = p
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    # local freestream accelerations
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate aerodynamic loads
    N, A, M = quasisteady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return SVector(L, M)
end

# convenience function for defining the coupling function parameters
function quasisteady_section_set_parameters!(p; U, rho, c)
    p[1] = U
    p[2] = rho
    p[3] = c
    return p
end
    
# convenience function for separating the coupling function parameters
quasisteady_section_separate_parameters(p) = (U = p[1], rho = p[2], c = p[3])

