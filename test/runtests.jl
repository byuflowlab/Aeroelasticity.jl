using AerostructuralDynamics
using GXBeam
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

@testset "Lifting Line" begin

    model = LiftingLine{3}(Wagner())

    dλ = rand(number_of_states(model))
    λ = rand(number_of_states(model))
    d = rand(number_of_inputs(model))
    p = rand(number_of_parameters(model))
    t = rand()

    fr = (λ) -> get_rates(model, λ, d, p, t)
    fin = (d) -> get_rates(model, λ, d, p, t)

    J = AerostructuralDynamics.get_state_jacobian(model, λ, d, p, t)

    Jy = AerostructuralDynamics.get_input_jacobian(model, λ, d, p, t)

    # test core jacobian functions
    @test isapprox(ForwardDiff.jacobian(fr, λ), J)
    @test isapprox(ForwardDiff.jacobian(fin, d), Array(Jy))
end

@testset "Geometrically Exact Beam Theory" begin

    # number of beam elements
    nelem = 2

    # points which define beam elements
    x = range(0, 60, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # index of endpoints of each beam element
    start = 1:nelem
    stop = 2:nelem+1

    # stiffness matrix for each beam element
    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        nelem)

    # mass matrix for each beam element
    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         nelem)

    # create assembly of interconnected nonlinear beams
    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

    # set prescribed conditions
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            )

    # set distributed loads
    distributed = Dict(i => DistributedLoads(assembly, i) for i = 1:nelem)

    # create system
    system = System(assembly, keys(prescribed), false)

    # create structural model
    model = GEBT(system, assembly, prescribed, distributed)

    # extract system constants and  pointers
    force_scaling = system.force_scaling
    mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam

    # set origin, linear velocity, and angular velocity
    x0 = zeros(3)
    v0 = zeros(3)
    ω0 = zeros(3)

    # get element properties and indices at which distributed loads are applied
    elements = assembly.elements
    element_indices = keys(distributed)

    # sample state rates, states, loads, parameters, and time
    dq = rand(number_of_states(model))
    q = rand(number_of_states(model))
    r = rand(number_of_inputs(model))
    p = rand(number_of_parameters(model))
    t = rand()

    # function for LHS, state rate input
    fl = (dq) -> AerostructuralDynamics.gxbeam_lhs!(similar(dq), q, dq, assembly, prescribed, distributed,
        force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2,
        icol_pt, icol_beam, x0, v0, ω0)

    # function for RHS, state input
    fr = (q) -> AerostructuralDynamics.get_rates(model::GEBT, q, r, p, t)

    # function for RHS, load input
    fin = (r) -> AerostructuralDynamics.get_rates(model::GEBT, q, r, p, t)

    # storage for mass and jacobian matrices
    M = similar(q, length(q), length(q))
    J = similar(q, length(q), length(q))

    # calculate mass matrix
    AerostructuralDynamics.gxbeam_mass_matrix!(M, q, assembly, force_scaling,
        mass_scaling, irow_pt, irow_beam, irow_beam1, irow_beam2, icol_pt,
        icol_beam)

    # calculate jacobian with respect to states
    AerostructuralDynamics.gxbeam_state_jacobian!(J, q, assembly, prescribed,
        distributed, force_scaling, mass_scaling, irow_pt, irow_beam,
        irow_beam1, irow_beam2, icol_pt, icol_beam, x0, v0, ω0)

    # calculate jacobian with respect to loads
    Jy = AerostructuralDynamics.gxbeam_input_jacobian(q, r, elements, start,
        stop, force_scaling, mass_scaling, irow_pt, element_indices)

    # test mass matrix
    @test all(isapprox.(ForwardDiff.jacobian(fl, dq), M, atol=1e-10))

    # test state jacobian
    @test all(isapprox.(ForwardDiff.jacobian(fr, q), J, atol=1e-10))

    # test input jacobian
    @test all(isapprox.(ForwardDiff.jacobian(fin, r), Array(Jy), atol=1e-10))

end
