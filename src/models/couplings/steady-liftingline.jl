
"""
    SteadyLiftingLine

Coupling model which allows a steady thin airfoil theory model (see [`Steady`](@ref)) to
be used more conveniently with the [`LiftingLine`](@ref) model. See [`LiftingLineSection`](@ref).
"""
struct SteadyLiftingLine
    steady::Steady
    liftingline::LiftingLineSection
end

default_coupling(steady::Steady, liftingline::LiftingLineSection) = SteadyLiftingLine(steady, liftingline)

function (steady_liftingline::SteadyLiftingLine)(dx, x, p, t)
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = x[2]
    # extract parameters
    a, b, a0, α0, cd0, cm0 = p[1]
    ρ, c = p[2]
    # freestream velocity components
    u, v, ω = liftingline_section_velocities(vx, vz, -ωy)
    # calculate aerodynamic loads
    N, A, M = steady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return nothing, (f, m)
end