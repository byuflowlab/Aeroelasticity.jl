
"""
    QuasiSteadyLiftingLine

Coupling model which allows a quasisteady thin airfoil theory model (see [`Steady`](@ref)) 
to be used more conveniently with the [`LiftingLine`](@ref) model. See [`LiftingLineSection`](@ref).
"""
struct QuasiSteadyLiftingLine
    quasisteady::QuasiSteady
    liftingline::LiftingLineSection
end

default_coupling(quasisteady::QuasiSteady, liftingline::LiftingLineSection) = QuasiSteadyLiftingLine(quasisteady, liftingline)

function (quasisteady_liftingline::QuasiSteadyLiftingLine)(dx, x, p, t)
    # rate variables
    dvx, dvy, dvz, dωx, dωy, dωz = dx[2]
    # extract state variables
    vx, vy, vz, ωx, ωy, ωz = x[2]
    # extract parameters
    a, b, a0, α0, cd0, cm0 = p[1]
    ρ, c = p[2]
    # freestream velocity components
    u, v, ω = liftingline_section_velocities(vx, vz, ωy)
    udot, vdot, ωdot = liftingline_section_accelerations(dvx, dvz, dωy)
    # calculate aerodynamic loads
    N, A, M = quasisteady_loads(a, b, ρ, c, a0, α0, cd0, cm0, u, v, ω, vdot, ωdot)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return nothing, (f, m)
end
