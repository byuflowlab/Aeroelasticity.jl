# --- Coupling Model Creation --- #

"""
    WagnerSection

Coupling model for coupling an unsteady aerodynamic model based on Wagner's function (see
[`Wagner`](@ref)) and a two-degree of freedom typical section model (see [`Section`](@ref)).
This coupling introduces the freestream velocity ``U_\\infty``, air density
``\\rho_\\infty``, and the Prandtl-Glauert compressibility factor ``\\beta`` as additional 
parameters.

The parameters for the resulting coupled model (as defined by the parameter function)
defaults to the parameters for each model concatenated into a single vector.
"""
struct WagnerSection{TF}
    wagner::Wagner{TF}
    section::Section
end

default_coupling(wagner::Wagner, section::Section) = WagnerSection(wagner, section)

function (wagner_section::WagnerSection)(dx, x, p, t)
    # extract constants
    @unpack C1, C2 = wagner_section.wagner
    # extract rate variables
    dλ1, dλ2 = dx[1]
    dh, dθ, dhdot, dθdot = dx[2]
    # extract state variables
    λ1, λ2 = x[1]
    h, θ, hdot, θdot = x[2]
    # extract parameters
    a, b, a0, α0, cd0, cm0 = p[1]
    kh, kθ, m, Sθ, Iθ = p[2]
    U, rho, beta = p[3]
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    N, A, M = wagner_loads(a, b, a0, α0, cd0, cm0, rho, beta, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # lift is approximately normal force
    L = N
    # return inputs
    return SVector(u, v, ω), SVector(L, M)
end

# default parameter function
default_parameter_function(::WagnerSection) = (p, t) -> (view(p, 1:6), view(p, 7:11), view(p, 12:14))