# --- Coupling Model Creation --- #

"""
    QuasiSteadySection

Coupling model for coupling a quasi-steady aerodynamic model based on thin airfoil theory
(see [`QuasiSteady`](@ref)) and a two-degree of freedom typical section model
(see [`Section()`]).  This model introduces the freestream velocity ``U``, air density
``\\rho``, and the Prandtl-Glauert compressibility factor ``\\beta`` as additional 
parameters.

The parameters for the resulting coupled model (as defined by the parameter function)
defaults to the parameters for each model concatenated into a single vector.
"""
struct QuasiSteadySection
    quasisteady::QuasiSteady
    section::Section
end

default_coupling(quasisteady::QuasiSteady, section::Section) = QuasiSteadySection(quasisteady, section)

function (::QuasiSteadySection)(dx, x, p, t)
    # extract state variables
    dh, dθ, dhdot, dθdot = dx[2]
    # extract state variables
    h, θ, hdot, θdot = x[2]
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, cd0, cm0 = p[1]
    kh, kθ, m, Sθ, Iθ = p[2]
    U, rho, beta = p[3]
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    # local freestream accelerations
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate aerodynamic loads
    N, A, M = quasisteady_loads(a, b, a0, α0, cd0, cm0, rho, beta, u, v, ω, vdot, ωdot)
    # lift is approximately equal to the normal force
    L = N
    # return inputs
    return nothing, SVector(L, M)
end

# default parameter function
default_parameter_function(::QuasiSteadySection) = (p, t) -> (view(p, 1:6), view(p, 7:11), view(p, 12:14))
