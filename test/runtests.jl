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

@testset "Peter's Finite State" begin

    model = PetersFiniteState{4}()

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

@testset "Peters Finite State" begin

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

@testet "Peters Finite State - Hodges" begin

    using AerostructuralDynamics, LinearAlgebra

    # number of aerodynamic states
    N = 6

    # dimensionless parameters
    a = -1/5 # center of mass normalized location
    e = -1/10 # reference point normalized location
    μ = 20 # = m/(ρ*pi*b^2)
    r2 = 6/25 # = Ip/(m*b^2)
    σ = 2/5 # = ωh/ωθ
    V = range(0, 3, length=100)

    # dimensionalized parameters
    b = 0.5
    ρ = 1
    m = ρ*pi*b^2
    ωθ = 1
    ωh = ωθ*σ
    kh = m*ωh^2
    kθ = m*ωθ^2
    xθ = e - a
    Ip = r2*m*b^2
    U = V*b*ωθ

    # aerodynamic model
    aero = PetersFiniteState{N}()

    # structural model
    stru = TypicalSection()

    # combined models
    models = (aero, stru)

    # --- Test Aerodynamic Model --- #

    U0 = 0
    u_aero = zeros(N)
    y_aero = zeros(N)
    p_aero = [a, b, U0, ρ]
    t = 0

    dλ = get_rates(aero, u_aero, y_aero, p_aero, t)
    Jλ = get_state_jacobian(aero, u_aero, y_aero, p_aero, t)
    Mλ = get_mass_matrix(aero, u_aero, y_aero, p_aero, t)

    # --- Test Structural Model --- #

    u_stru = zeros(4)
    y_stru = zeros(2)
    p_stru = [a, b, kh, kθ, m, xθ, Ip]
    t = 0

    dq = get_rates(stru, u_stru, y_stru, p_stru, t)
    Jq = get_state_jacobian(stru, u_stru, y_stru, p_stru, t)
    Mq = get_mass_matrix(stru, u_stru, y_stru, p_stru, t)

    @test Mq == [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
    @test Jq == [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]

    # test zero velocity eigenvalues
    u_aero = zeros(N)
    u_stru = zeros(4)
    u = vcat(u_aero, u_stru)

    # parameters
    p_aero = [a, b, U[i], ρ]
    p_stru = [a, b, kh, kθ, m, xθ, Ip]
    p = vcat(p_aero, p_stru)

    # time
    t = 0.0

    # calculate inputs
    y = get_inputs(models, u, p, t)




    get_mass_matrix(stru, )
end
