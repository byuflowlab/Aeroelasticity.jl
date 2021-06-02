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
    a, b, kh, kθ, m, xθ, Ip = p

    fl = (dq) -> AerostructuralDynamics.section_lhs(b, m, xθ, Ip, dq...)
    fr = (q) -> AerostructuralDynamics.section_rhs(a, b, kh, kθ, q..., L, M)
    fin = (r) -> AerostructuralDynamics.section_rhs(a, b, kh, kθ, q..., r...)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fl, dq),
        AerostructuralDynamics.section_mass_matrix(b, m, xθ, Ip))
    @test isapprox(ForwardDiff.jacobian(fr, q),
        AerostructuralDynamics.section_state_jacobian(kh, kθ))
    @test isapprox(ForwardDiff.jacobian(fin, r),
        AerostructuralDynamics.section_input_jacobian(a, b))

    # TODO: Test interface functions
end

@testset "QuasiSteady" begin

    dq = rand(4)
    q = rand(4)
    p = rand(4)

    dh, dθ, dhdot, dθdot = dq
    h, θ, hdot, θdot = q
    a, b, U, ρ = p

    fout0_q = (q) -> AerostructuralDynamics.zero_order_loads(b, U, ρ, q[2])
    fout1_q = (q) -> fout0_q(q) .+ AerostructuralDynamics.first_order_loads(a, b,
        U, ρ, q[3], q[4])
    fout2_q = (q) -> fout1_q(q) .+ AerostructuralDynamics.second_order_loads(a, b, U, ρ, q[4], 0, 0)
    fout2_dq = (dq) -> fout1_q(q) .+ AerostructuralDynamics.second_order_loads(a, b, U, ρ, 0, dq[3], dq[4])

    @test isapprox(ForwardDiff.jacobian(fout0_q, q),
        AerostructuralDynamics.quasisteady0_jacobian(b, U, ρ))
    @test isapprox(ForwardDiff.jacobian(fout1_q, q),
        AerostructuralDynamics.quasisteady1_jacobian(a, b, U, ρ))
    @test isapprox(ForwardDiff.jacobian(fout2_q, q),
        AerostructuralDynamics.quasisteady2_jacobian(a, b, U, ρ))
    @test isapprox(-ForwardDiff.jacobian(fout2_dq, dq),
        AerostructuralDynamics.quasisteady2_mass_matrix(a, b, U, ρ))

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
    θ, hdot, θdot = d
    a, b, U, ρ = p

    C1 = model.C1
    C2 = model.C2
    ε1 = model.eps1
    ε2 = model.eps2

    fr = (λ) -> AerostructuralDynamics.wagner_rates(a, b, U, C1, C2, ε1, ε2, θ, hdot, θdot, λ...)
    fin = (d) -> AerostructuralDynamics.wagner_rates(a, b, U, C1, C2, ε1, ε2, d..., λ1, λ2)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fr, λ),
        AerostructuralDynamics.wagner_state_jacobian(a, b, U, ε1, ε2))
    @test isapprox(ForwardDiff.jacobian(fin, d),
        AerostructuralDynamics.wagner_input_jacobian(a, b, U, C1, C2, ε1, ε2))

    h = rand()
    θ = rand()
    hdot = rand()
    θdot = rand()
    hddot = rand()
    θddot = rand()

    fout_λ = (λ) -> AerostructuralDynamics.wagner_loads(a, b, U, ρ, C1, C2, θ,
        hdot, θdot, hddot, θddot, λ...)
    fout_h = (h) -> AerostructuralDynamics.wagner_loads(a, b, U, ρ, C1, C2, θ,
        hdot, θdot, hddot, θddot, λ1, λ2)
    fout_θ = (θ) -> AerostructuralDynamics.wagner_loads(a, b, U, ρ, C1, C2, θ,
        hdot, θdot, hddot, θddot, λ1, λ2)
    fout_hdot = (hdot) -> AerostructuralDynamics.wagner_loads(a, b, U, ρ, C1, C2, θ,
        hdot, θdot, hddot, θddot, λ1, λ2)
    fout_θdot = (θdot) -> AerostructuralDynamics.wagner_loads(a, b, U, ρ, C1, C2, θ,
        hdot, θdot, hddot, θddot, λ1, λ2)
    fout_hddot = (hddot) -> AerostructuralDynamics.wagner_loads(a, b, U, ρ, C1, C2, θ,
        hdot, θdot, hddot, θddot, λ1, λ2)
    fout_θddot = (θddot) -> AerostructuralDynamics.wagner_loads(a, b, U, ρ, C1, C2, θ,
        hdot, θdot, hddot, θddot, λ1, λ2)

    @test isapprox(ForwardDiff.jacobian(fout_λ, λ),
        AerostructuralDynamics.wagner_loads_λ(b, U, ρ))
    @test isapprox(ForwardDiff.derivative(fout_h, h),
        AerostructuralDynamics.wagner_loads_h())
    @test isapprox(ForwardDiff.derivative(fout_θ, θ),
        AerostructuralDynamics.wagner_loads_θ(b, U, ρ, C1, C2))
    @test isapprox(ForwardDiff.derivative(fout_hdot, hdot),
        AerostructuralDynamics.wagner_loads_hdot(b, U, ρ, C1, C2))
    @test isapprox(ForwardDiff.derivative(fout_θdot, θdot),
        AerostructuralDynamics.wagner_loads_θdot(a, b, U, ρ, C1, C2))
    @test isapprox(ForwardDiff.derivative(fout_hddot, hddot),
        AerostructuralDynamics.wagner_loads_hddot(b, ρ))
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

    θdot, hddot, θddot = d
    a, b, U, ρ = p

    fl = (dλ) -> AerostructuralDynamics.peters_lhs(model.A, dλ)
    fr = (λ) -> AerostructuralDynamics.peters_rhs(a, b, U, model.c, θdot, hddot, θddot, λ)
    fin = (d) -> AerostructuralDynamics.peters_rhs(a, b, U, model.c, d..., λ)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fl, dλ),
        AerostructuralDynamics.peters_mass_matrix(model.A))
    @test isapprox(ForwardDiff.jacobian(fr, λ),
        AerostructuralDynamics.peters_state_jacobian(b, U, model.c))
    @test isapprox(ForwardDiff.jacobian(fin, d),
        AerostructuralDynamics.peters_input_jacobian(a, b, U, model.c))

    h = rand()
    θ = rand()
    hdot = rand()
    θdot = rand()
    hddot = rand()
    θddot = rand()

    fout_λ = (λ) -> AerostructuralDynamics.peters_loads(a, b, U, ρ, model.b, θ,
        hdot, θdot, hddot, θddot, λ)
    fout_h = (h) -> AerostructuralDynamics.peters_loads(a, b, U, ρ, model.b, θ,
        hdot, θdot, hddot, θddot, λ)
    fout_θ = (θ) -> AerostructuralDynamics.peters_loads(a, b, U, ρ, model.b, θ,
        hdot, θdot, hddot, θddot, λ)
    fout_hdot = (hdot) -> AerostructuralDynamics.peters_loads(a, b, U, ρ, model.b, θ,
        hdot, θdot, hddot, θddot, λ)
    fout_θdot = (θdot) -> AerostructuralDynamics.peters_loads(a, b, U, ρ, model.b, θ,
        hdot, θdot, hddot, θddot, λ)
    fout_hddot = (hddot) -> AerostructuralDynamics.peters_loads(a, b, U, ρ, model.b, θ,
        hdot, θdot, hddot, θddot, λ)
    fout_θddot = (θddot) -> AerostructuralDynamics.peters_loads(a, b, U, ρ, model.b, θ,
        hdot, θdot, hddot, θddot, λ)

    @test isapprox(ForwardDiff.jacobian(fout_λ, λ),
        AerostructuralDynamics.peters_loads_λ(b, U, ρ, model.b))
    @test isapprox(ForwardDiff.derivative(fout_h, h),
        AerostructuralDynamics.peters_loads_h())
    @test isapprox(ForwardDiff.derivative(fout_θ, θ),
        AerostructuralDynamics.peters_loads_θ(b, U, ρ))
    @test isapprox(ForwardDiff.derivative(fout_hdot, hdot),
        AerostructuralDynamics.peters_loads_hdot(b, U, ρ))
    @test isapprox(ForwardDiff.derivative(fout_θdot, θdot),
        AerostructuralDynamics.peters_loads_θdot(a, b, U, ρ))
    @test isapprox(ForwardDiff.derivative(fout_hddot, hddot),
        AerostructuralDynamics.peters_loads_hddot(b, ρ))
    @test isapprox(ForwardDiff.derivative(fout_θddot, θddot),
        AerostructuralDynamics.peters_loads_θddot(a, b, ρ))

    # TODO: Test interface functions
end
