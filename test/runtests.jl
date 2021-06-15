using AerostructuralDynamics
using LinearAlgebra
using GXBeam
using ForwardDiff
using Test

function run_model_tests(model; atol = 1e-10, norm = (x)->norm(x, Inf))

    du = rand(number_of_states(model))
    u = rand(number_of_states(model))
    y = rand(number_of_inputs(model))
    p = rand(number_of_parameters(model))
    t = rand()

    # # left hand side, state rates
    # fdu = (du) -> get_mass_matrix_product(model, du, u, y, p, t)

    # right hand side, state variables
    fu = (u) -> get_rates(model, u, y, p, t)

    # right hand side, inputs
    fy = (y) -> get_rates(model, u, y, p, t)

    # residual, state rates
    fresid = (du) -> get_mass_matrix(model, u, y, p, t)*du - get_rates(model, u, y, p, t)

    # # mass matrix test
    # M = get_mass_matrix(model, u, y, p, t)
    # M_fd = ForwardDiff.jacobian(fdu, du)
    # @test isapprox(M, M_fd; atol, norm)

    # state jacobian test
    J = get_state_jacobian(model, u, y, p, t)
    J_fd = ForwardDiff.jacobian(fu, u)
    @test isapprox(J, J_fd; atol, norm)

    # input jacobian test
    Jy = get_input_jacobian(model, u, y, p, t)
    Jy_fd = ForwardDiff.jacobian(fy, y)
    @test isapprox(Array(Jy), Jy_fd; atol, norm)

    # # compatability test
    # du = qr(collect(M))\get_rates(model, u, y, p, t)
    # lhs = get_mass_matrix_product(model, du, u, y, p, t)
    # rhs = get_rates(models, u, y, p, t)
    # @test isapprox(lhs, rhs; atol, norm)

    return nothing
end

function run_model_tests(models...; atol = 1e-10, norm = (x)->norm(x, Inf))

    du = rand(number_of_states(models))
    u = rand(number_of_states(models))
    p = rand(number_of_parameters(models))
    t = rand()

    y = get_inputs(models, u, p, t)

    # left hand side, state rates
    fdu = (du) -> get_mass_matrix_product(models, du, u, get_inputs(models, u, p, t), p, t)

    # right hand side, state variables
    fu = (u) -> get_rates(models, u, get_inputs(models, u, p, t), p, t)

    # right hand side, inputs
    fy = (y) -> get_rates(models, u, y, p, t)

    # residual, state rates
    fresid = (du) -> get_mass_matrix(models, u,
        get_inputs(models, u, p, t), p, t)*du - get_rates(models, u,
        get_inputs(models, u, p, t), p, t)

    # # mass matrix test
    # M = get_mass_matrix(models, u, y, p, t)
    # M_fd = ForwardDiff.jacobian(fdu, du)
    # @test isapprox(M, M_fd; atol, norm)

    # state jacobian test
    J = get_state_jacobian(models, u, y, p, t)
    J_fd = ForwardDiff.jacobian(fu, u)
    @test isapprox(J, J_fd; atol, norm)

    # input jacobian test
    Jy = get_input_jacobian(models, u, y, p, t)
    Jy_fd = ForwardDiff.jacobian(fy, y)
    @test isapprox(Array(Jy), Jy_fd; atol, norm)

    # # compatability test
    # du = qr(collect(M))\get_rates(models, u, y, p, t)
    # lhs = get_mass_matrix_product(models, du, u, y, p, t)
    # rhs = get_rates(models, u, y, p, t)
    # @test isapprox(lhs, rhs; atol, norm)

    return nothing
end

@testset "Typical Section" begin
    # test on its own
    run_model_tests(TypicalSection())
end

@testset "QuasiSteady" begin
    # test on its own
    run_model_tests(QuasiSteady{0}())
    run_model_tests(QuasiSteady{1}())
    run_model_tests(QuasiSteady{2}())
    # test coupling with TypicalSection()
    run_model_tests(QuasiSteady{0}(), TypicalSection())
    run_model_tests(QuasiSteady{1}(), TypicalSection())
    run_model_tests(QuasiSteady{2}(), TypicalSection())
end

@testset "Wagner" begin
    # test on its own
    run_model_tests(Wagner())
    # test coupling with TypicalSection()
    run_model_tests(Wagner(), TypicalSection())
end

@testset "Peters' Finite State" begin
    # test on its own
    run_model_tests(Peters{4}())
    # test coupling with TypicalSection()
    run_model_tests(Peters{4}(), TypicalSection())
end

@testset "Geometrically Exact Beam Theory" begin

    # test on its own

    N = 4

    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    start = 1:N
    stop = 2:N+1

    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        N)

    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         N)

    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            )

    distributed = Dict(i => DistributedLoads(assembly, i) for i = 1:N)

    system = System(assembly, keys(prescribed), false)

    model = GEBT(system, assembly, prescribed, distributed)

    run_model_tests(model)
end

@testset "Lifting Line" begin

    # number of lifting line sections
    N = 4

    # test on its own
    run_model_tests(LiftingLine{N}(QuasiSteady{0}()))
    run_model_tests(LiftingLine{N}(QuasiSteady{1}()))
    run_model_tests(LiftingLine{N}(QuasiSteady{2}()))
    run_model_tests(LiftingLine{N}(Wagner()))
    run_model_tests(LiftingLine{N}(Peters{4}()))

    # test coupling with GXBeam

    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    start = 1:N
    stop = 2:N+1

    stiffness = fill(
        [2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8],
        N)

    mass = fill(
        [258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418],
         N)

    assembly = Assembly(points, start, stop; stiffness=stiffness, mass=mass)

    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            )

    distributed = Dict(i => DistributedLoads(assembly, i) for i = 1:N)

    system = System(assembly, keys(prescribed), false)

    structural_model = GEBT(system, assembly, prescribed, distributed)

    run_model_tests(LiftingLine{4}(QuasiSteady{0}()), structural_model)
    run_model_tests(LiftingLine{4}(QuasiSteady{1}()), structural_model)
    run_model_tests(LiftingLine{4}(QuasiSteady{2}()), structural_model)
    run_model_tests(LiftingLine{4}(Wagner()), structural_model)
    run_model_tests(LiftingLine{4}(Peters{4}()), structural_model)
end
