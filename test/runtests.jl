using AerostructuralDynamics
using LinearAlgebra
using GXBeam
using ForwardDiff
using Test

const AD = AerostructuralDynamics

function run_model_jacobian_tests(model;
    dx = rand(number_of_states(model)),
    x = rand(number_of_states(model)),
    y = rand(number_of_inputs(model)),
    p = rand(number_of_parameters(model)),
    t = rand(),
    atol = sqrt(eps()),
    norm = (x) -> LinearAlgebra.norm(x, Inf))

    # rate jacobian test
    f = (dx) -> get_residual(model, dx, x, y, p, t)
    J = get_rate_jacobian(model, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, dx)
    @test isapprox(J, J_fd; atol, norm)

    # state jacobian test
    f = (x) -> get_residual(model, dx, x, y, p, t)
    J = get_state_jacobian(model, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, x)
    @test isapprox(J, J_fd; atol, norm)

    # input jacobian test
    f = (y) -> get_residual(model, dx, x, y, p, t)
    J = get_input_jacobian(model, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, y)
    @test isapprox(J, J_fd; atol, norm)

    # parameter jacobian test
    f = (p) -> get_residual(model, dx, x, y, p, t)
    J = get_parameter_jacobian(model, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, p)
    @test isapprox(J, J_fd; atol, norm)

    # time derivative test
    f = (t) -> get_residual(model, dx, x, y, p, t)
    dT = get_time_gradient(model, dx, x, y, p, t)
    dT_fd = ForwardDiff.derivative(f, t)
    @test isapprox(dT, dT_fd; atol, norm)

    return nothing
end

function run_coupling_jacobian_tests(models...;
    dx = rand(number_of_states(models)),
    x = rand(number_of_states(models)),
    p = rand(number_of_parameters(models)),
    t = rand(),
    atol = sqrt(eps()),
    norm = (x)->norm(x, Inf))

    # rate jacobian test
    f = (dx) -> get_coupling_inputs(models, dx, x, p, t)
    J = AerostructuralDynamics.get_coupling_rate_jacobian(models, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, dx)
    @test isapprox(J, J_fd; atol, norm)

    # state jacobian test
    f = (x) -> get_coupling_inputs(models, dx, x, p, t)
    J = AerostructuralDynamics.get_coupling_state_jacobian(models, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, x)
    @test isapprox(J, J_fd; atol, norm)

    # parameter jacobian test
    f = (p) -> get_coupling_inputs(models, dx, x, p, t)
    J = AerostructuralDynamics.get_coupling_parameter_jacobian(models, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, p)
    @test isapprox(J, J_fd; atol, norm)

    # time derivative test
    f = (t) -> get_coupling_inputs(models, dx, x, p, t)
    dT = AerostructuralDynamics.get_coupling_time_gradient(models, dx, x, p, t)
    dT_fd = ForwardDiff.derivative(f, t)
    @test isapprox(dT, dT_fd; atol, norm)

    return nothing
end

function run_model_io_tests(model; atol = sqrt(eps()), norm = (x)->norm(x, Inf))

    # test separating and recombining state variables
    x = rand(number_of_states(model))
    states = separate_states(model, x)
    new_x = get_states(model; states...)
    @test isapprox(new_x, x; atol, norm)

    # test separating and recombining inputs
    y = rand(number_of_inputs(model))
    inputs = separate_inputs(model, y)
    new_y = get_inputs(model; inputs...)
    @test isapprox(new_y, y; atol, norm)

    # test separating and recombining parameters
    p = rand(number_of_parameters(model))
    parameters = separate_parameters(model, p)
    new_p = get_parameters(model; parameters...)
    @test isapprox(new_p, p; atol, norm)

    return nothing
end

function run_coupling_io_tests(models...; atol = sqrt(eps()), norm = (x)->norm(x, Inf))

    # test separating and recombining additional parameters
    padd = rand(AerostructuralDynamics.number_of_additional_parameters(models...))
    parameters = AerostructuralDynamics.separate_additional_parameters(models..., padd)
    new_padd = get_additional_parameters(models; parameters...)
    @test isapprox(new_padd, padd; atol, norm)

    return nothing
end

# --- Two-Dimensional Models --- #

# Aerodynamic Models (2D)

@testset "QuasiSteady" begin
    model = QuasiSteady{0}()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

    model = QuasiSteady{1}()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

    model = QuasiSteady{2}()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Wagner" begin
    model = Wagner()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Peters' Finite State" begin
    model = Peters{4}()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

# Structural Models (2D)

@testset "Typical Section" begin
    model = TypicalSection()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

# Control Surface Models (2D)

@testset "SimpleFlap" begin
    model = SimpleFlap()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

# Controller Models (2D)

# Coupled Models (2D)
@testset "QuasiSteady + TypicalSection" begin
    models = (QuasiSteady{0}(), TypicalSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{1}(), TypicalSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{2}(), TypicalSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Wagner + TypicalSection" begin
    models = (Wagner(), TypicalSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Peters' Finite State + TypicalSection" begin
    models = (Peters{4}(), TypicalSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "QuasiSteady + TypicalSection + SimpleFlap" begin
    models = (QuasiSteady{0}(), TypicalSection(), SimpleFlap())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{1}(), TypicalSection(), SimpleFlap())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{2}(), TypicalSection(), SimpleFlap())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Wagner + TypicalSection + SimpleFlap" begin
    models = (Wagner(), TypicalSection(), SimpleFlap())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Peters + TypicalSection + SimpleFlap" begin
    models = (Peters{4}(), TypicalSection(), SimpleFlap())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "QuasiSteady + LiftingLineSection" begin
    models = (QuasiSteady{0}(), AerostructuralDynamics.LiftingLineSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{1}(), AerostructuralDynamics.LiftingLineSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{2}(), AerostructuralDynamics.LiftingLineSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Wagner + LiftingLineSection" begin
    models = (Wagner(), AerostructuralDynamics.LiftingLineSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Peters' Finite State + LiftingLineSection" begin
    models = (Peters{4}(), AerostructuralDynamics.LiftingLineSection())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "QuasiSteady + LiftingLineSection + SimpleFlap + LiftingLineSectionControl" begin
    models = (QuasiSteady{0}(), AerostructuralDynamics.LiftingLineSection(),
        SimpleFlap(), AerostructuralDynamics.LiftingLineSectionControl())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{1}(), AerostructuralDynamics.LiftingLineSection(),
        SimpleFlap(), AerostructuralDynamics.LiftingLineSectionControl())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (QuasiSteady{2}(), AerostructuralDynamics.LiftingLineSection(),
        SimpleFlap(), AerostructuralDynamics.LiftingLineSectionControl())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Wagner + LiftingLineSection + SimpleFlap + LiftingLineSectionControl" begin
    models = (Wagner(), AerostructuralDynamics.LiftingLineSection(),
        SimpleFlap(), AerostructuralDynamics.LiftingLineSectionControl())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Peters + LiftingLineSection + SimpleFlap + LiftingLineSectionControl" begin
    models = (Peters{4}(), AerostructuralDynamics.LiftingLineSection(),
        SimpleFlap(), AerostructuralDynamics.LiftingLineSectionControl())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

# --- Three Dimensional Models --- #

# Aerodynamic Models (3D)

@testset "Lifting Line" begin

    model = LiftingLine{2}(QuasiSteady{0}())
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

    model = LiftingLine{2}(QuasiSteady{1}())
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

    model = LiftingLine{2}(QuasiSteady{2}())
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

    model = LiftingLine{2}(Wagner())
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

    model = LiftingLine{2}(Peters{4}())
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

end

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

    # run model tests
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

# Dynamics Models (3D)

@testset "RigidBody" begin

    model = RigidBody()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

end

# Control Surface Models (3D)

@testset "LiftingLineFlaps" begin

    N = 2 # number of lifting line sections
    NG = 1 # number of control inputs
    model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

end

# Controller Models (3D)

@testset "Trim" begin

    model = Trim()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)

end

# Coupled Models

@testset "Lifting Line + RigidBody" begin

    # number of lifting line sections
    N = 2

    models = (LiftingLine{N}(QuasiSteady{0}()), RigidBody())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (LiftingLine{N}(QuasiSteady{1}()), RigidBody())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (LiftingLine{N}(QuasiSteady{2}()), RigidBody())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (LiftingLine{N}(Wagner()), RigidBody())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    models = (LiftingLine{N}(Peters{4}()), RigidBody())
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

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

    # define models
    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    structural_model = GEBT(assembly, prescribed)
    models = (aerodynamic_model, structural_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    structural_model = GEBT(assembly, prescribed)
    models = (aerodynamic_model, structural_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    structural_model = GEBT(assembly, prescribed)
    models = (aerodynamic_model, structural_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    structural_model = GEBT(assembly, prescribed)
    models = (aerodynamic_model, structural_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    structural_model = GEBT(assembly, prescribed)
    models = (aerodynamic_model, structural_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
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

    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    models = (aerodynamic_model, structural_model, dynamics_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    models = (aerodynamic_model, structural_model, dynamics_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    models = (aerodynamic_model, structural_model, dynamics_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    models = (aerodynamic_model, structural_model, dynamics_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    models = (aerodynamic_model, structural_model, dynamics_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

end

@testset "Lifting Line + RigidBody + Lifting Line Flaps" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 1

    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

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

    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

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

    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
end

@testset "Lifting Line + RigidBody + Lifting Line Flaps + Trim" begin

    # number of lifting line sections
    N = 2

    # number of control surface inputs
    NG = 6

    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
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

    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    structural_model = GEBT(assembly, prescribed)
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)
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

    aerodynamic_model = LiftingLine{N}(QuasiSteady{0}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{1}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(QuasiSteady{2}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Wagner())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

    aerodynamic_model = LiftingLine{N}(Peters{4}())
    structural_model = GEBT(assembly, prescribed)
    dynamics_model = RigidBody()
    surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
    control_model = Trim()
    models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
    run_coupling_jacobian_tests(models...)
    run_coupling_io_tests(models...)

end
