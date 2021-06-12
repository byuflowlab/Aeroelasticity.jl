using AerostructuralDynamics
using ForwardDiff
using Test

@testset "Typical Section" begin

    model = TypicalSection()

    dq = rand(number_of_states(model))
    q = rand(number_of_states(model))
    r = rand(number_of_inputs(model))
    p = rand(number_of_parameters(model))
    t = rand()

    dh, dθ, dhdot, dθdot = dq
    h, θ, hdot, θdot = q
    L, M = r
    kh, kθ, m, Sθ, Iθ = p

    fl = (dq) -> AerostructuralDynamics.section_lhs(m, Sθ, Iθ, dq...)
    fr = (q) -> AerostructuralDynamics.section_rhs(kh, kθ, q..., L, M)
    fin = (r) -> AerostructuralDynamics.section_rhs(kh, kθ, q..., r...)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fl, dq),
        AerostructuralDynamics.section_mass_matrix(m, Sθ, Iθ))
    @test isapprox(ForwardDiff.jacobian(fr, q),
        AerostructuralDynamics.section_state_jacobian(kh, kθ))
    @test isapprox(ForwardDiff.jacobian(fin, r),
        AerostructuralDynamics.section_input_jacobian())

    # TODO: Test interface functions
end

@testset "QuasiSteady" begin

    dq = rand(4)
    q = rand(4)
    p = rand(5)

    dh, dθ, dhdot, dθdot = dq
    h, θ, hdot, θdot = q
    a, b, ρ, a0, α0 = p
    u = rand()

    v(u, θ, hdot) = -u*θ - hdot
    fout0_q = (q) -> AerostructuralDynamics.quasisteady0_loads(a, b, ρ, a0,
        α0, u, v(u, q[2], 0))
    fout1_q = (q) -> AerostructuralDynamics.quasisteady1_loads(a, b, ρ, a0,
        α0, u, v(u, q[2], q[3]), q[4])
    fout2_q = (q) -> AerostructuralDynamics.quasisteady2_loads(a, b, ρ, a0,
        α0, u, v(u, q[2], q[3]), -dhdot, q[4], dθdot)
    fout2_dq = (dq) -> AerostructuralDynamics.quasisteady2_loads(a, b, ρ, a0,
        α0, u, v(u, θ, hdot), -dq[3], θdot, dq[4])

    @test isapprox(ForwardDiff.jacobian(fout0_q, q),
        AerostructuralDynamics.quasisteady0_jacobian(a, b, ρ, a0, u))
    @test isapprox(ForwardDiff.jacobian(fout1_q, q),
        AerostructuralDynamics.quasisteady1_jacobian(a, b, ρ, a0, u))
    @test isapprox(ForwardDiff.jacobian(fout2_q, q),
        AerostructuralDynamics.quasisteady2_jacobian(a, b, ρ, a0, u))
    @test isapprox(-ForwardDiff.jacobian(fout2_dq, dq),
        AerostructuralDynamics.quasisteady2_mass_matrix(a, b, ρ))

end

@testset "Wagner" begin

    model = Wagner()

    dλ = rand(number_of_states(model))
    λ = rand(number_of_states(model))
    d = rand(number_of_inputs(model))
    p = rand(number_of_parameters(model))
    t = rand()

    dλ1, dλ2 = dλ
    λ1, λ2 = λ
    u, v, θdot = d
    a, b, ρ, a0, α0 = p

    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2

    fr = (λ) -> AerostructuralDynamics.wagner_rates(a, b, α0, C1, C2, ε1, ε2,
        u, v, θdot, λ...)
    fin = (d) -> AerostructuralDynamics.wagner_rates(a, b, α0, C1, C2, ε1, ε2,
        d..., λ1, λ2)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fr, λ),
        AerostructuralDynamics.wagner_state_jacobian(a, b, ε1, ε2, u))
    @test isapprox(ForwardDiff.jacobian(fin, d),
        AerostructuralDynamics.wagner_input_jacobian(a, b, α0, C1, C2, ε1, ε2,
        u, v, θdot, λ1, λ2))

    h = rand()
    θ = rand()
    hdot = rand()
    θdot = rand()
    hddot = rand()
    θddot = rand()

    fout_λ = (λ) -> AerostructuralDynamics.wagner_loads(a, b, ρ, a0, α0, C1,
        C2, u, -u*θ - hdot, -hddot, θdot, θddot, λ...)
    fout_h = (h) -> AerostructuralDynamics.wagner_loads(a, b, ρ, a0, α0, C1,
        C2, u, -u*θ - hdot, -hddot, θdot, θddot, λ1, λ2)
    fout_θ = (θ) -> AerostructuralDynamics.wagner_loads(a, b, ρ, a0, α0, C1,
        C2, u, -u*θ - hdot, -hddot, θdot, θddot, λ1, λ2)
    fout_hdot = (hdot) -> AerostructuralDynamics.wagner_loads(a, b, ρ, a0,
        α0, C1, C2, u, -u*θ - hdot, -hddot, θdot, θddot, λ1, λ2)
    fout_θdot = (θdot) -> AerostructuralDynamics.wagner_loads(a, b, ρ, a0,
        α0, C1, C2, u, -u*θ - hdot, -hddot, θdot, θddot, λ1, λ2)
    fout_hddot = (hddot) -> AerostructuralDynamics.wagner_loads(a, b, ρ, a0,
        α0, C1, C2, u, -u*θ - hdot, -hddot, θdot, θddot, λ1, λ2)
    fout_θddot = (θddot) -> AerostructuralDynamics.wagner_loads(a, b, ρ,
        a0, α0, C1, C2, u, -u*θ - hdot, -hddot, θdot, θddot, λ1, λ2)

    @test isapprox(ForwardDiff.jacobian(fout_λ, λ),
        AerostructuralDynamics.wagner_loads_λ(a, b, ρ, a0, u))
    @test isapprox(ForwardDiff.derivative(fout_h, h),
        AerostructuralDynamics.wagner_loads_h())
    @test isapprox(ForwardDiff.derivative(fout_θ, θ),
        AerostructuralDynamics.wagner_loads_θ(a, b, ρ, a0, C1, C2, u))
    @test isapprox(ForwardDiff.derivative(fout_hdot, hdot),
        AerostructuralDynamics.wagner_loads_hdot(a, b, ρ, a0, C1, C2, u))
    @test isapprox(ForwardDiff.derivative(fout_θdot, θdot),
        AerostructuralDynamics.wagner_loads_θdot(a, b, ρ, a0, C1, C2, u))
    @test isapprox(ForwardDiff.derivative(fout_hddot, hddot),
        AerostructuralDynamics.wagner_loads_hddot(a, b, ρ))
    @test isapprox(ForwardDiff.derivative(fout_θddot, θddot),
        AerostructuralDynamics.wagner_loads_θddot(a, b, ρ))

    # TODO: Test interface functions
end

@testset "Peters' Finite State" begin

    model = Peters{4}()

    dλ = rand(number_of_states(model))
    λ = rand(number_of_states(model))
    d = rand(number_of_inputs(model))
    p = rand(number_of_parameters(model))
    t = rand()

    u, vdot, θdot, θddot = d
    a, b, ρ, a0, α0 = p

    fl = (dλ) -> AerostructuralDynamics.peters_lhs(model.A, dλ)
    fr = (λ) -> AerostructuralDynamics.peters_rhs(a, b, model.c, u, vdot, θdot, θddot, λ)
    fin = (d) -> AerostructuralDynamics.peters_rhs(a, b, model.c, d..., λ)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fl, dλ),
        AerostructuralDynamics.peters_mass_matrix(model.A))
    @test isapprox(ForwardDiff.jacobian(fr, λ),
        AerostructuralDynamics.peters_state_jacobian(b, model.c, u))
    @test isapprox(ForwardDiff.jacobian(fin, d),
        AerostructuralDynamics.peters_input_jacobian(a, b, model.c, u, θdot, λ))

    h = rand()
    θ = rand()
    hdot = rand()
    θdot = rand()
    hddot = rand()
    θddot = rand()

    fout_λ = (λ) -> AerostructuralDynamics.peters_loads(a, b, ρ, a0, α0,
        model.b, u, -u*θ - hdot, -hddot, θdot, θddot, λ)
    fout_h = (h) -> AerostructuralDynamics.peters_loads(a, b, ρ, a0, α0,
        model.b, u, -u*θ - hdot, -hddot, θdot, θddot, λ)
    fout_θ = (θ) -> AerostructuralDynamics.peters_loads(a, b, ρ, a0, α0,
        model.b, u, -u*θ - hdot, -hddot, θdot, θddot, λ)
    fout_hdot = (hdot) -> AerostructuralDynamics.peters_loads(a, b, ρ, a0,
        α0, model.b, u, -u*θ - hdot, -hddot, θdot, θddot, λ)
    fout_θdot = (θdot) -> AerostructuralDynamics.peters_loads(a, b, ρ, a0,
        α0, model.b, u, -u*θ - hdot, -hddot, θdot, θddot, λ)
    fout_hddot = (hddot) -> AerostructuralDynamics.peters_loads(a, b, ρ, a0,
        α0, model.b, u, -u*θ - hdot, -hddot, θdot, θddot, λ)
    fout_θddot = (θddot) -> AerostructuralDynamics.peters_loads(a, b, ρ, a0,
        α0, model.b, u, -u*θ - hdot, -hddot, θdot, θddot, λ)

    @test isapprox(ForwardDiff.jacobian(fout_λ, λ),
        AerostructuralDynamics.peters_loads_λ(a, b, ρ, a0, model.b, u))
    @test isapprox(ForwardDiff.derivative(fout_h, h),
        AerostructuralDynamics.peters_loads_h())
    @test isapprox(ForwardDiff.derivative(fout_θ, θ),
        AerostructuralDynamics.peters_loads_θ(a, b, ρ, a0, u))
    @test isapprox(ForwardDiff.derivative(fout_hdot, hdot),
        AerostructuralDynamics.peters_loads_hdot(a, b, ρ, a0, u))
    @test isapprox(ForwardDiff.derivative(fout_θdot, θdot),
        AerostructuralDynamics.peters_loads_θdot(a, b, ρ, a0, u))
    @test isapprox(ForwardDiff.derivative(fout_hddot, hddot),
        AerostructuralDynamics.peters_loads_hddot(a, b, ρ))
    @test isapprox(ForwardDiff.derivative(fout_θddot, θddot),
        AerostructuralDynamics.peters_loads_θddot(a, b, ρ))

    # TODO: Test interface functions
end
