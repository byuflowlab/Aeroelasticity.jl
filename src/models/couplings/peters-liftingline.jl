"""
    PetersLiftingLine{N,TF,TV,TA}

Coupling model which allows Peters' finite state theory (see [`Peters`](@ref)) to be 
used more conveniently with the [`LiftingLine`](@ref) model. See [`LiftingLineSection`](@ref).
"""
struct PetersLiftingLine{N,TF,TV,TA}
    peters::Peters{N,TF,TV,TA}
    liftingline::LiftingLineSection
end

default_coupling(peters::Peters, liftingline::LiftingLineSection) = PetersLiftingLine(peters, liftingline)

function (peters_liftingline::PetersLiftingLine)(dx, x, p, t)
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
    ρ, c = p[2]
    # local freestream velocity components
    u, v, ω = liftingline_section_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_section_accelerations(dvx, dvz, dωy)
    # calculate loads
    N, A, M = peters_loads(a, b, ρ, c, a0, α0, cd0, cm0, bbar, u, v, ω, vdot, ωdot, λ)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs for each model
    return SVector(u, ω, vdot, ωdot), (f, m)
end