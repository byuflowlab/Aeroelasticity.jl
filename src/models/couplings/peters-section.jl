
"""
    PetersSection{N}

Coupling model for coupling Peters' finite state theory (see [`Peters`](@ref)) with a
typical section model (see [`Section`](@ref)).  This model introduces the freestream
velocity ``U``, air density ``\\rho``, and the Prandtl-Glauert compressibility factor 
``\\beta`` as additional parameters.

The parameters for the resulting coupled model (as defined by the parameter function)
defaults to the parameters for each model concatenated into a single vector.
"""
struct PetersSection{N,TF,TV,TA}
    peters::Peters{N,TF,TV,TA}
    section::Section
end

default_coupling(peters::Peters, section::Section) = PetersSection(peters, section)

function (peters_section::PetersSection)(dx, x, p, t)
    # extract constants
    bbar = peters_section.peters.b
    # extract rate variables
    dλ = dx[1]
    dh, dθ, dhdot, dθdot = dx[2]
    # extract state variables
    λ = x[1]
    h, θ, hdot, θdot = x[2]
    # extract parameters
    a, b, a0, α0, cd0, cm0 = p[1]
    kh, kθ, m, Sθ, Iθ = p[2]
    U, rho, beta = p[3]
    # local freestream velocity components
    u, v, ω = section_velocities(U, θ, hdot, θdot)
    udot, vdot, ωdot = section_accelerations(dhdot, dθdot)
    # calculate loads
    N, A, M = peters_loads(a, b, a0, α0, cd0, cm0, rho, beta, bbar, u, v, ω, vdot, ωdot, λ)
    # lift is approximately equal to the normal force
    L = N
    # return inputs for each model
    return SVector(u, ω, vdot, ωdot), SVector(L, M)
end

# default parameter function
default_parameter_function(::PetersSection) = (p, t) -> (view(p, 1:6), view(p, 7:11), view(p, 12:14))
