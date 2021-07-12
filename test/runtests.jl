using AerostructuralDynamics
using LinearAlgebra
using GXBeam
using ForwardDiff
using Test

const AD = AerostructuralDynamics

function has_mass_matrix(model)
    matrix_type = AerostructuralDynamics.mass_matrix_type(typeof(model))
    empty_matrix = AerostructuralDynamics.isempty(matrix_type)
    identity_matrix = AerostructuralDynamics.isidentity(matrix_type)
    return !(empty_matrix || identity_matrix)
end

function has_input_mass_matrix(models...)
    input_matrix_type = AerostructuralDynamics.mass_matrix_type(typeof.(models)...)
    zero_matrix = AerostructuralDynamics.iszero(input_matrix_type)
    return !zero_matrix
end

function run_model_tests(model; atol = 1e-10, norm = (x)->norm(x, Inf),
    test_mass_matrix = has_mass_matrix(model))

    # sample arguments
    du = rand(number_of_states(model))
    u = rand(number_of_states(model))
    y = rand(number_of_inputs(model))
    p = rand(number_of_parameters(model))
    t = rand()

    # state jacobian test
    fu = (u) -> get_rates(model, u, y, p, t)
    J = get_state_jacobian(model, u, y, p, t)
    J_fd = ForwardDiff.jacobian(fu, u)
    @test isapprox(J, J_fd; atol, norm)

    # input jacobian test
    fy = (y) -> get_rates(model, u, y, p, t)
    Jy = get_input_jacobian(model, u, y, p, t)
    Jy_fd = ForwardDiff.jacobian(fy, y)
    @test isapprox(Array(Jy), Jy_fd; atol, norm)

    # additional tests if left hand side is defined
    if test_mass_matrix

        # mass matrix test
        fdu = (du) -> AD.get_lhs(model, du, u, y, p, t)
        M = get_mass_matrix(model, u, y, p, t)
        M_fd = ForwardDiff.jacobian(fdu, similar(u))
        @test isapprox(M, M_fd; atol, norm)

        # compatability test
        du = qr(collect(M), Val(true))\get_rates(model, u, y, p, t)
        lhs = AD.get_lhs(model, du, u, y, p, t)
        @test isapprox(M*du, lhs; atol, norm)

    end

    return nothing
end

function run_coupling_tests(models...; atol = 1e-10, norm = (x)->norm(x, Inf),
    test_mass_matrix = has_input_mass_matrix(models...))

    # sample arguments
    du = rand(number_of_states(models))
    u = rand(number_of_states(models))
    p = rand(number_of_parameters(models))
    t = rand()

    # jacobian test
    fu = (u) -> get_inputs(models, u, p, t)
    Jy = AD.get_input_state_jacobian(models, u, p, t)
    Jy_fd = ForwardDiff.jacobian(fu, u)
    @test isapprox(Jy, Jy_fd; atol, norm)

    # additional test if state rate contribution is defined
    if test_mass_matrix

        # mass matrix test
        fdu = (du) -> -AD.get_inputs_from_state_rates(models..., du, u, p, t)
        My = AD.get_input_mass_matrix(models, u, p, t)
        My_fd = ForwardDiff.jacobian(fdu, du)
        @test isapprox(My, My_fd; atol, norm)

    end

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
    run_coupling_tests(QuasiSteady{0}(), TypicalSection())
    run_coupling_tests(QuasiSteady{1}(), TypicalSection())
    run_coupling_tests(QuasiSteady{2}(), TypicalSection())
    # test coupling with LiftingLineSection()
    run_coupling_tests(QuasiSteady{0}(), AD.LiftingLineSection())
    run_coupling_tests(QuasiSteady{1}(), AD.LiftingLineSection())
    run_coupling_tests(QuasiSteady{2}(), AD.LiftingLineSection())
end

@testset "Wagner" begin
    # test on its own
    run_model_tests(Wagner())
    # test coupling with TypicalSection()
    run_coupling_tests(Wagner(), TypicalSection())
    # test coupling with LiftingLineSection()
    run_coupling_tests(Wagner(), AD.LiftingLineSection())
end

@testset "Peters' Finite State" begin
    # test on its own
    run_model_tests(Peters{4}())
    # test coupling with TypicalSection()
    run_coupling_tests(Peters{4}(), TypicalSection())
    # test coupling with LiftingLineSection()
    run_coupling_tests(Peters{4}(), AD.LiftingLineSection())
end

@testset "Rigid Body" begin
    run_model_tests(RigidBody())
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

    # test coupling with RigidBody
    run_coupling_tests(LiftingLine{4}(QuasiSteady{0}()), RigidBody())
    run_coupling_tests(LiftingLine{4}(QuasiSteady{1}()), RigidBody())
    run_coupling_tests(LiftingLine{4}(QuasiSteady{2}()), RigidBody())
    run_coupling_tests(LiftingLine{4}(Wagner()), RigidBody())
    run_coupling_tests(LiftingLine{4}(Peters{4}()), RigidBody())

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

    run_coupling_tests(LiftingLine{4}(QuasiSteady{0}()), structural_model)
    run_coupling_tests(LiftingLine{4}(QuasiSteady{1}()), structural_model)
    run_coupling_tests(LiftingLine{4}(QuasiSteady{2}()), structural_model)
    run_coupling_tests(LiftingLine{4}(Wagner()), structural_model)
    run_coupling_tests(LiftingLine{4}(Peters{4}()), structural_model)

    run_coupling_tests(LiftingLine{4}(QuasiSteady{0}()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{4}(QuasiSteady{1}()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{4}(QuasiSteady{2}()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{4}(Wagner()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{4}(Peters{4}()), structural_model, RigidBody())
end
