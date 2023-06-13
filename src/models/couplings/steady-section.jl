# --- Coupling Model Creation --- #

"""
    SteadySection

Coupling model for coupling a steady aerodynamic model based on thin airfoil theory
(see [`Steady`](@ref)) and a two-degree of freedom typical section model
(see [`Section`](@ref)).  This model introduces the freestream velocity ``U_\\infty``, air
density ``\\rho_\\infty``, and the Prandtl-Glauert compressibility factor ``\\beta`` as 
additional parameters.

The parameters for the resulting coupled model (as defined by the parameter function)
defaults to the parameters for each model concatenated into a single vector.
"""
struct SteadySection
    steady::Steady
    section::Section
end

default_coupling(steady::Steady, section::Section) = SteadySection(steady, section)

function (::SteadySection)(dx, x, p, t)
    # extract state variables
    h, θ, hdot, θdot = x[2]
    # extract aerodynamic, structural, and aerostructural parameters
    a, b, a0, α0, cd0, cm0 = p[1]
    kh, kθ, m, Sθ, Iθ = p[2]
    U, rho, beta = p[3]
    # local freestream velocity components
    u, v = section_steady_velocities(U, θ)
    # calculate aerodynamic loads
    N, A, M = steady_loads(a, b, a0, α0, cd0, cm0, rho, beta, u, v)
    # lift is approximately equal to the normal force
    L = N
    # return inputs for each model
    return nothing, SVector(L, M)
end

# default parameter function
default_parameter_function(::SteadySection) = (p, t) -> (view(p, 1:6), view(p, 7:11), view(p, 12:14))

# --- Internal Methods --- #

# local freestream velocities
function section_steady_velocities(U, θ)
    u = U
    v = U*θ
    return SVector(u, v)
end