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

function has_coupling_mass_matrix(models...)
    coupling_matrix_type = AerostructuralDynamics.coupling_mass_matrix_type(typeof.(models)...)
    zero_matrix = AerostructuralDynamics.iszero(coupling_matrix_type)
    return !zero_matrix
end

function run_model_tests(model;
    dx = rand(number_of_states(model)),
    x = rand(number_of_states(model)),
    y = rand(number_of_inputs(model)),
    p = rand(number_of_parameters(model)),
    t = rand(),
    atol = sqrt(eps()),
    norm = (x)->norm(x, Inf),
    test_mass_matrix = has_mass_matrix(model),
    test_convenience_functions = true)

    # state jacobian test
    fu = (x) -> get_rates(model, x, y, p, t)
    J = get_state_jacobian(model, x, y, p, t)
    J_fd = ForwardDiff.jacobian(fu, x)
    @test isapprox(J, J_fd; atol, norm)

    # input jacobian test
    fy = (y) -> get_rates(model, x, y, p, t)
    Jy = AD.get_input_jacobian(model, x, y, p, t)
    Jy_fd = ForwardDiff.jacobian(fy, y)
    @test isapprox(Array(Jy), Jy_fd; atol, norm)

    # additional tests if left hand side is defined
    if test_mass_matrix

        # mass matrix test
        fdu = (dx) -> AD.get_lhs(model, dx, x, y, p, t)
        M = get_mass_matrix(model, x, y, p, t)
        M_fd = ForwardDiff.jacobian(fdu, similar(x))
        @test isapprox(M, M_fd; atol, norm)

        # compatability test
        dx = qr(collect(M), Val(true))\get_rates(model, x, y, p, t)
        lhs = AD.get_lhs(model, dx, x, y, p, t)
        @test isapprox(M*dx, lhs; atol, norm)

    end

    # convenience function tests
    if test_convenience_functions
        states = separate_states(model, x)
        inputs = separate_inputs(model, y)
        parameters = separate_parameters(model, p)
        new_x = get_states(model; states...)
        new_y = get_inputs(model; inputs...)
        new_p = get_parameters(model; parameters...)
        @test isapprox(new_x, x)
        @test isapprox(new_y, y)
        @test isapprox(new_p, p)
    end

    return nothing
end

function run_coupling_tests(models...;
    dx = rand(number_of_states(models)),
    x = rand(number_of_states(models)),
    p = rand(number_of_parameters(models)),
    t = rand(),
    atol = sqrt(eps()),
    norm = (x)->norm(x, Inf),
    test_mass_matrix = has_coupling_mass_matrix(models...),
    test_convenience_functions = true)

    # jacobian test
    fu = (x) -> get_inputs(models, x, p, t)
    Jy = AD.get_coupling_state_jacobian(models, x, p, t)
    Jy_fd = ForwardDiff.jacobian(fu, x)
    @test isapprox(Jy, Jy_fd; atol, norm)

    # additional test if state rate contribution is defined
    if test_mass_matrix
        fdu = (dx) -> -AD.get_inputs_using_state_rates(models..., dx, x, p, t)
        My = AD.get_coupling_mass_matrix(models, x, p, t)
        My_fd = ForwardDiff.jacobian(fdu, dx)
        @test isapprox(My, My_fd; atol, norm)
    end

    if test_convenience_functions
        padd = rand(AD.number_of_additional_parameters(models...))
        parameters = AD.separate_additional_parameters(models..., padd)
        new_padd = get_additional_parameters(models; parameters...)
        @test isapprox(new_padd, padd)
    end

    return nothing
end

# # --- Two-Dimensional Models --- #
#
# # Aerodynamic Models (2D)
#
# @testset "QuasiSteady" begin
#     # test on its own
#     run_model_tests(QuasiSteady{0}())
#     run_model_tests(QuasiSteady{1}())
#     run_model_tests(QuasiSteady{2}())
# end
#
# @testset "Wagner" begin
#     # test on its own
#     run_model_tests(Wagner())
# end
#
# @testset "Peters' Finite State" begin
#     # test on its own
#     run_model_tests(Peters{4}())
# end
#
# # Structural Models (2D)
#
# @testset "Typical Section" begin
#     # test on its own
#     run_model_tests(TypicalSection())
# end
#
# # Control Surface Models (2D)
#
# @testset "LinearFlap" begin
#     run_model_tests(LinearFlap())
# end
#
# # Controller Models (2D)
#
# # Coupled Models (2D)
# @testset "QuasiSteady + TypicalSection" begin
#     # test coupling with TypicalSection()
#     run_coupling_tests(QuasiSteady{0}(), TypicalSection())
#     run_coupling_tests(QuasiSteady{1}(), TypicalSection())
#     run_coupling_tests(QuasiSteady{2}(), TypicalSection())
# end
#
# @testset "Wagner + TypicalSection" begin
#     # test coupling with TypicalSection()
#     run_coupling_tests(Wagner(), TypicalSection())
# end
#
# @testset "Peters' Finite State + TypicalSection" begin
#     # test coupling with TypicalSection()
#     run_coupling_tests(Peters{4}(), TypicalSection())
# end
#
# @testset "QuasiSteady + TypicalSection + LinearFlap" begin
#     run_coupling_tests(QuasiSteady{0}(), TypicalSection(), LinearFlap())
#     run_coupling_tests(QuasiSteady{1}(), TypicalSection(), LinearFlap())
#     run_coupling_tests(QuasiSteady{2}(), TypicalSection(), LinearFlap())
# end
#
# @testset "Wagner + TypicalSection + LinearFlap" begin
#     run_coupling_tests(Wagner(), TypicalSection(), LinearFlap())
# end
#
# @testset "Peters + TypicalSection + LinearFlap" begin
#     run_coupling_tests(Peters{4}(), TypicalSection(), LinearFlap())
# end
#
# @testset "QuasiSteady + LiftingLineSection" begin
#     # test coupling with LiftingLineSection()
#     run_coupling_tests(QuasiSteady{0}(), AD.LiftingLineSection())
#     run_coupling_tests(QuasiSteady{1}(), AD.LiftingLineSection())
#     run_coupling_tests(QuasiSteady{2}(), AD.LiftingLineSection())
# end
#
# @testset "Wagner + LiftingLineSection" begin
#     # test coupling with LiftingLineSection()
#     run_coupling_tests(Wagner(), AD.LiftingLineSection())
# end
#
# @testset "Peters' Finite State + LiftingLineSection" begin
#     # test coupling with LiftingLineSection()
#     run_coupling_tests(Peters{4}(), AD.LiftingLineSection())
# end
#
# @testset "QuasiSteady + LiftingLineSection + LinearFlap + LiftingLineSectionControl" begin
#     run_coupling_tests(QuasiSteady{0}(), AD.LiftingLineSection(), LinearFlap(), AD.LiftingLineSectionControl())
#     run_coupling_tests(QuasiSteady{1}(), AD.LiftingLineSection(), LinearFlap(), AD.LiftingLineSectionControl())
#     run_coupling_tests(QuasiSteady{2}(), AD.LiftingLineSection(), LinearFlap(), AD.LiftingLineSectionControl())
# end
#
# @testset "Wagner + LiftingLineSection + LinearFlap + LiftingLineSectionControl" begin
#     run_coupling_tests(Wagner(), AD.LiftingLineSection(), LinearFlap(), AD.LiftingLineSectionControl())
# end
#
# @testset "Peters + LiftingLineSection + LinearFlap + LiftingLineSectionControl" begin
#     run_coupling_tests(Peters{4}(), AD.LiftingLineSection(), LinearFlap(), AD.LiftingLineSectionControl())
# end
#
# # --- Three Dimensional Models --- #
#
# # Aerodynamic Models (3D)
#
# @testset "Lifting Line" begin
#
#     # number of lifting line sections
#     N = 2
#
#     # test on its own
#     run_model_tests(LiftingLine{N}(QuasiSteady{0}()))
#     run_model_tests(LiftingLine{N}(QuasiSteady{1}()))
#     run_model_tests(LiftingLine{N}(QuasiSteady{2}()))
#     run_model_tests(LiftingLine{N}(Wagner()))
#     run_model_tests(LiftingLine{N}(Peters{4}()))
# end

# Structural Models (3D)

@testset "Geometrically Exact Beam Theory" begin

    # number of elements
    nelem = 4

    # point locations
    x = range(0, 60, length=nelem+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # element connectivity
    start = 1:nelem
    stop = 2:nelem+1

    # element orientation
    e1 = [1, 0, 0]
    e2 = [0, 1, 0]
    e3 = [0, 0, 1]
    frame = hcat(e1, e2, e3)

    # element stiffness properties
    stiffness = [
         2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
        ]
    compliance = inv(stiffness)
    compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

    # element inertial properties
    mass = [
         258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418
         ]

    # construct assembly
    assembly = Assembly(points, start, stop;
        frames = fill(frame, nelem),
        stiffness = fill(stiffness, nelem),
        mass = fill(mass, nelem))

    # define boundary conditions and applied loads
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            nelem+1 => PrescribedConditions(Fz=1e5)
            )

    # define model
    model = GEBT(assembly, prescribed)

    # # use provided parameters (instead of random parameters)
    # p = default_parameters(model, assembly)

    # # use provided inputs (instead of random inputs)
    # y = default_inputs(model, assembly; prescribed)

    run_model_tests(model)
end

# Dynamics Models (3D)

@testset "RigidBody" begin
    run_model_tests(RigidBody())
end

# Control Surface Models (3D)

@testset "LiftingLineFlaps" begin

    # number of lifting line sections
    N = 2

    # number of control inputs
    NG = 1

    # test on its own
    run_model_tests(LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))
end

# Controller Models (3D)

@testset "Trim" begin
    run_model_tests(Trim())
end

# Coupled Models

@testset "Lifting Line + RigidBody" begin

    # number of lifting line sections
    N = 2

    # test coupling with RigidBody
    run_coupling_tests(LiftingLine{N}(QuasiSteady{0}()), RigidBody())
    run_coupling_tests(LiftingLine{N}(QuasiSteady{1}()), RigidBody())
    run_coupling_tests(LiftingLine{N}(QuasiSteady{2}()), RigidBody())
    run_coupling_tests(LiftingLine{N}(Wagner()), RigidBody())
    run_coupling_tests(LiftingLine{N}(Peters{4}()), RigidBody())
end

@testset "Lifting Line + GXBeam" begin

    # number of lifting line sections
    N = 2

    # point locations
    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # element connectivity
    start = 1:N
    stop = 2:N+1

    # element orientation
    e1 = [1, 0, 0]
    e2 = [0, 1, 0]
    e3 = [0, 0, 1]
    frame = hcat(e1, e2, e3)

    # element stiffness properties
    stiffness = [
         2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
        ]
    compliance = inv(stiffness)
    compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

    # element inertial properties
    mass = [
         258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418
         ]

    # construct assembly
    assembly = Assembly(points, start, stop;
        frames = fill(frame, N),
        stiffness = fill(stiffness, N),
        mass = fill(mass, N))

    # define boundary conditions and applied loads
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            N+1 => PrescribedConditions(Fz=1e5)
            )

    # define model
    structural_model = GEBT(assembly, prescribed)

    run_coupling_tests(LiftingLine{N}(QuasiSteady{0}()), structural_model)
    run_coupling_tests(LiftingLine{N}(QuasiSteady{1}()), structural_model)
    run_coupling_tests(LiftingLine{N}(QuasiSteady{2}()), structural_model)
    run_coupling_tests(LiftingLine{N}(Wagner()), structural_model)
    run_coupling_tests(LiftingLine{N}(Peters{4}()), structural_model)
end

@testset "Lifting Line + GXBeam + RigidBody" begin

    # number of lifting line sections
    N = 2

    # point locations
    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # element connectivity
    start = 1:N
    stop = 2:N+1

    # element orientation
    e1 = [1, 0, 0]
    e2 = [0, 1, 0]
    e3 = [0, 0, 1]
    frame = hcat(e1, e2, e3)

    # element stiffness properties
    stiffness = [
         2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
        ]
    compliance = inv(stiffness)
    compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

    # element inertial properties
    mass = [
         258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418
         ]

    # construct assembly
    assembly = Assembly(points, start, stop;
        frames = fill(frame, N),
        stiffness = fill(stiffness, N),
        mass = fill(mass, N))

    # define boundary conditions and applied loads
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            N+1 => PrescribedConditions(Fz=1e5)
            )

    # define model
    structural_model = GEBT(assembly, prescribed)

    run_coupling_tests(LiftingLine{N}(QuasiSteady{0}()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{N}(QuasiSteady{1}()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{N}(QuasiSteady{2}()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{N}(Wagner()), structural_model, RigidBody())
    run_coupling_tests(LiftingLine{N}(Peters{4}()), structural_model, RigidBody())
end

@testset "Lifting Line + RigidBody + Lifting Line Flaps" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 1

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{0}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{1}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{2}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(Wagner()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(Peters{4}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))
end

@testset "Lifting Line + GXBeam + Lifting Line Flaps" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 1

    # point locations
    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # element connectivity
    start = 1:N
    stop = 2:N+1

    # element orientation
    e1 = [1, 0, 0]
    e2 = [0, 1, 0]
    e3 = [0, 0, 1]
    frame = hcat(e1, e2, e3)

    # element stiffness properties
    stiffness = [
         2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
        ]
    compliance = inv(stiffness)
    compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

    # element inertial properties
    mass = [
         258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418
         ]

    # construct assembly
    assembly = Assembly(points, start, stop;
        frames = fill(frame, N),
        stiffness = fill(stiffness, N),
        mass = fill(mass, N))

    # define boundary conditions and applied loads
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            N+1 => PrescribedConditions(Fz=1e5)
            )

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{0}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{1}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{2}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(Wagner()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(Peters{4}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))
end

@testset "Lifting Line + GXBeam + RigidBody + Lifting Line Flaps" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 1

    # point locations
    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # element connectivity
    start = 1:N
    stop = 2:N+1

    # element orientation
    e1 = [1, 0, 0]
    e2 = [0, 1, 0]
    e3 = [0, 0, 1]
    frame = hcat(e1, e2, e3)

    # element stiffness properties
    stiffness = [
         2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
        ]
    compliance = inv(stiffness)
    compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

    # element inertial properties
    mass = [
         258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418
         ]

    # construct assembly
    assembly = Assembly(points, start, stop;
        frames = fill(frame, N),
        stiffness = fill(stiffness, N),
        mass = fill(mass, N))

    # define boundary conditions and applied loads
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            N+1 => PrescribedConditions(Fz=1e5)
            )

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{0}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{1}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{2}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(Wagner()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))

    run_coupling_tests(
        LiftingLine{N}(Peters{4}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)))
end

@testset "Lifting Line + RigidBody + Lifting Line Flaps + Trim" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 6

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{0}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{1}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{2}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(Wagner()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(Peters{4}()),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())
end

@testset "Lifting Line + GXBeam + Lifting Line Flaps + Trim" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 6

    # point locations
    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # element connectivity
    start = 1:N
    stop = 2:N+1

    # element orientation
    e1 = [1, 0, 0]
    e2 = [0, 1, 0]
    e3 = [0, 0, 1]
    frame = hcat(e1, e2, e3)

    # element stiffness properties
    stiffness = [
         2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
        ]
    compliance = inv(stiffness)
    compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

    # element inertial properties
    mass = [
         258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418
         ]

    # construct assembly
    assembly = Assembly(points, start, stop;
        frames = fill(frame, N),
        stiffness = fill(stiffness, N),
        mass = fill(mass, N))

    # define boundary conditions and applied loads
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            N+1 => PrescribedConditions(Fz=1e5)
            )

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{0}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{1}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{2}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(Wagner()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(Peters{4}()),
        GEBT(assembly, prescribed),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())
end

@testset "Lifting Line + GXBeam + RigidBody + Lifting Line Flaps + Trim" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 6

    # point locations
    x = range(0, 60, length=N+1)
    y = zero(x)
    z = zero(x)
    points = [[x[i],y[i],z[i]] for i = 1:length(x)]

    # element connectivity
    start = 1:N
    stop = 2:N+1

    # element orientation
    e1 = [1, 0, 0]
    e2 = [0, 1, 0]
    e3 = [0, 0, 1]
    frame = hcat(e1, e2, e3)

    # element stiffness properties
    stiffness = [
         2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
         1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
         6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
        -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
        -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
        -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
        ]
    compliance = inv(stiffness)
    compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

    # element inertial properties
    mass = [
         258.053      0.0        0.0      0.0      7.07839  -71.6871
           0.0      258.053      0.0     -7.07839  0.0        0.0
           0.0        0.0      258.053   71.6871   0.0        0.0
           0.0       -7.07839   71.6871  48.59     0.0        0.0
           7.07839    0.0        0.0      0.0      2.172      0.0
         -71.6871     0.0        0.0      0.0      0.0       46.418
         ]

    # construct assembly
    assembly = Assembly(points, start, stop;
        frames = fill(frame, N),
        stiffness = fill(stiffness, N),
        mass = fill(mass, N))

    # define boundary conditions and applied loads
    prescribed = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            N+1 => PrescribedConditions(Fz=1e5)
            )

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{0}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{1}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(QuasiSteady{2}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(Wagner()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())

    run_coupling_tests(
        LiftingLine{N}(Peters{4}()),
        GEBT(assembly, prescribed),
        RigidBody(),
        LiftingLineFlaps{N}(LinearFlap(), ntuple(i->ones(N), NG)),
        Trim())
end
