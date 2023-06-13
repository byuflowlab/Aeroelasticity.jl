
"""
    WagnerLiftingLine{TF}

Coupling model which allows the Wagner unsteady aerodynamic model (see [`Wagner`](@ref)) to
be used more conveniently with the [`LiftingLine`](@ref) model. See [`LiftingLineSection`](@ref).
"""
struct WagnerLiftingLine{TF}
    wagner::Wagner{TF}
    liftingline::LiftingLineSection
end

default_coupling(wagner::Wagner, liftingline::LiftingLineSection) = WagnerLiftingLine(wagner, liftingline)

function (wagner_liftingline::WagnerLiftingLine)(dx, x, p, t)
    # extract constants
    @unpack C1, C2 = wagner_section.wagner
    # extract rate variables
    dλ1, dλ2 = dx[1]
    dvx, dvy, dvz, dωx, dωy, dωz = dx[2]
    # extract state variables
    λ1, λ2 = x[1]
    vx, vy, vz, ωx, ωy, ωz = x[2]
    # extract parameters
    a, b, a0, α0, cd0, cm0 = p[1]
    rho, beta = p[2]
    # freestream velocity components
    u, v, ω = liftingline_section_velocities(vx, vz, -ωy)
    udot, vdot, ωdot = liftingline_section_accelerations(dvx, dvz, -dωy)
    # calculate aerodynamic loads
    N, A, M = wagner_loads(a, b, a0, α0, cd0, cm0, rho, beta, C1, C2, u, v, ω, vdot, ωdot, λ1, λ2)
    # forces and moments per unit span
    f = SVector(A, 0, N)
    m = SVector(0, M, 0)
    # return inputs
    return SVector(u, v, ω), (f, m)
end