using AerostructuralDynamics
using ForwardDiff
using Test

# import functions used for testing
@testset "Typical Section" begin

    du = rand(4)
    u = rand(4)
    p = rand(7)
    r = rand(2)
    t = rand()

    dh, dθ, dhdot, dθdot = du
    h, θ, hdot, θdot = u
    a, b, kh, kθ, m, xθ, Ip = p
    L, M = r

    fl = (du) -> AerostructuralDynamics.section_lhs(b, m, xθ, Ip, du...)
    fr = (u) -> AerostructuralDynamics.section_rhs(a, b, kh, kθ, u..., L, M)
    fload = (r) -> AerostructuralDynamics.section_rhs(a, b, kh, kθ, u..., r...)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fl, du),
        AerostructuralDynamics.section_mass_matrix(b, m, xθ, Ip))
    @test isapprox(ForwardDiff.jacobian(fr, u),
        AerostructuralDynamics.section_state_jacobian(kh, kθ))
    @test isapprox(ForwardDiff.jacobian(fload, r),
        AerostructuralDynamics.section_load_jacobian(a, b))

    # TODO: Test interface functions

end

@testset "Peters Finite State + Typical Section" begin

    N = 4

    aero = PetersFiniteState{N}()
    stru = TypicalSection()

    dq = rand(4)
    dλ = rand(N)
    q = rand(4)
    λ = rand(N)
    p = rand(9)
    t = rand()

    dh, dθ, dhdot, dθdot = dq
    h, θ, hdot, θdot = q
    a, b, kh, kθ, m, xθ, Ip, U, ρ = p
    Abar, bbar, cbar = aero.A, aero.b, aero.c

    fstru_l = (dq) -> AerostructuralDynamics.peters_lhs(a, b, Abar, cbar, dq[3], dq[4], dλ)
    faero_l = (dλ) -> AerostructuralDynamics.peters_lhs(a, b, Abar, cbar, dhdot, dθdot, dλ)
    fstru_r = (q) -> AerostructuralDynamics.peters_rhs(b, U, cbar, q[4], λ)
    faero_r = (λ) -> AerostructuralDynamics.peters_rhs(b, U, cbar, θdot, λ)
    gstru_l = (q) -> -AerostructuralDynamics.peters_rate_loads(a, b, ρ, q[3], q[4])
    gaero_l = (λ) -> -AerostructuralDynamics.peters_rate_loads(a, b, ρ, dhdot, dθdot)
    gstru_r = (q) -> AerostructuralDynamics.peters_state_loads(a, b, U, ρ, bbar, q[2], q[3], q[4], λ)
    gaero_r = (λ) -> AerostructuralDynamics.peters_state_loads(a, b, U, ρ, bbar, θ, hdot, θdot, λ)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fstru_l, dq),
        AerostructuralDynamics.peters_stru_mass_matrix(a, b, cbar))
    @test isapprox(ForwardDiff.jacobian(faero_l, dλ),
        AerostructuralDynamics.peters_aero_mass_matrix(Abar))
    @test isapprox(ForwardDiff.jacobian(fstru_r, q),
        AerostructuralDynamics.peters_stru_jacobian(U, cbar))
    @test isapprox(ForwardDiff.jacobian(faero_r, λ),
        AerostructuralDynamics.peters_aero_jacobian(b, U, cbar))

    # test core jacobian functions for loads
    @test isapprox(ForwardDiff.jacobian(gstru_l, dq),
        AerostructuralDynamics.peters_load_stru_mass_matrix(a, b, ρ))
    @test isapprox(ForwardDiff.jacobian(gaero_l, dλ),
        AerostructuralDynamics.peters_load_aero_mass_matrix(cbar))
    @test isapprox(ForwardDiff.jacobian(gstru_r, q),
        AerostructuralDynamics.peters_load_stru_jacobian(a, b, U, ρ))
    @test isapprox(ForwardDiff.jacobian(gaero_r, λ),
        AerostructuralDynamics.peters_load_aero_jacobian(b, ρ, bbar))

    # TODO: Test interface functions

end
