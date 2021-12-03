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

function run_coupling_jacobian_tests(coupling;
    dx = rand(number_of_states(coupling)),
    x = rand(number_of_states(coupling)),
    p = rand(number_of_parameters(coupling)),
    t = rand(),
    atol = sqrt(eps()),
    norm = (x)->norm(x, Inf))

    # rate jacobian test
    f = (dx) -> get_coupling_inputs(coupling, dx, x, p, t)
    J = AerostructuralDynamics.get_coupling_rate_jacobian(coupling, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, dx)
    @test isapprox(J, J_fd; atol, norm)

    # state jacobian test
    f = (x) -> get_coupling_inputs(coupling, dx, x, p, t)
    J = AerostructuralDynamics.get_coupling_state_jacobian(coupling, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, x)
    @test isapprox(J, J_fd; atol, norm)

    # parameter jacobian test
    f = (p) -> get_coupling_inputs(coupling, dx, x, p, t)
    J = AerostructuralDynamics.get_coupling_parameter_jacobian(coupling, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, p)
    @test isapprox(J, J_fd; atol, norm)

    # time derivative test
    f = (t) -> get_coupling_inputs(coupling, dx, x, p, t)
    dT = AerostructuralDynamics.get_coupling_time_gradient(coupling, dx, x, p, t)
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

function run_coupling_io_tests(coupling; atol = sqrt(eps()), norm = (x)->norm(x, Inf))

    # test separating and recombining additional parameters
    pc = rand(AerostructuralDynamics.number_of_additional_parameters(coupling))
    parameters = AerostructuralDynamics.separate_parameters(coupling, pc)
    new_pc = get_parameters(coupling; parameters...)
    @test isapprox(new_pc, pc; atol, norm)

    return nothing
end

@testset "Steady" begin
    model = steady_model()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "QuasiSteady" begin
    model = quasisteady_model()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Wagner" begin
    model = wagner_model()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Peters' Finite State" begin
    model = peters_model(4)
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Typical Section" begin
    model = typical_section_model()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Lifting Line Section" begin
    model = liftingline_section_model()
    run_model_io_tests(model)
end

@testset "Steady + TypicalSection" begin
    coupled_model = steady_section_model()
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "QuasiSteady + TypicalSection" begin
    coupled_model = quasisteady_section_model()
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Wagner + TypicalSection" begin
    coupled_model = wagner_section_model()
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Peters' Finite State + TypicalSection" begin
    coupled_model = peters_section_model(4)
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Steady + Lifting Line Section" begin
    coupled_model = steady_liftingline_model()
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "QuasiSteady + Lifting Line Section" begin
    coupled_model = quasisteady_liftingline_model()
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Wagner + Lifting Line Section" begin
    coupled_model = wagner_liftingline_model()
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Peters' Finite State + Lifting Line Section" begin
    coupled_model = peters_liftingline_model(4)
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Steady Lifting Line" begin
    submodels = fill(steady_model(), 2)
    model = liftingline_model(submodels)
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Quasisteady Lifting Line" begin
    submodels = fill(quasisteady_model(), 2)
    model = liftingline_model(submodels)
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Wagner Lifting Line" begin
    submodels = fill(wagner_model(), 2)
    model = liftingline_model(submodels)
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Peters' Lifting Line" begin
    submodels = fill(peters_model(4), 2)
    model = liftingline_model(submodels)
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

# @testset "Geometrically Exact Beam Theory Initialization" begin
#     # number of elements
#     nelem = 2
#     # connectivity
#     start = 1:nelem
#     stop = 2:nelem+1
#     # boundary conditions
#     displacement = zeros(Bool, 6, nelem+1)
#     displacement[:,1] .= true # fixed at one end
#     # run model tests
#     model = gxbeam_model(start, stop, displacement)
#     run_model_jacobian_tests(model)
#     run_model_io_tests(model)
# end

@testset "Geometrically Exact Beam Theory" begin
    # number of elements
    nelem = 2
    # connectivity
    start = 1:nelem
    stop = 2:nelem+1
    # boundary conditions
    displacement = zeros(Bool, 6, nelem+1)
    displacement[:,1] .= true # fixed at one end
    # run model tests
    model = gxbeam_model(start, stop, displacement)
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Rigid Body" begin
    model = rigidbody_model()
    run_model_jacobian_tests(model)
    run_model_io_tests(model)
end

@testset "Steady Lifting Line + GXBeam" begin
    # number of elements
    nelem = 2
    # connectivity
    start = 1:nelem
    stop = 2:nelem+1
    # boundary conditions
    displacement = zeros(Bool, 6, nelem+1)
    displacement[:,1] .= true # fixed at one end
    # submodels
    lifting_elements = 1:nelem
    section_models = fill(steady_liftingline_model(), nelem)
    # define and test coupling
    coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
        section_models)
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Quasi-Steady Lifting Line + GXBeam" begin
    # number of elements
    nelem = 2
    # connectivity
    start = 1:nelem
    stop = 2:nelem+1
    # boundary conditions
    displacement = zeros(Bool, 6, nelem+1)
    displacement[:,1] .= true # fixed at one end
    # submodels
    lifting_elements = 1:nelem
    section_models = fill(quasisteady_liftingline_model(), nelem)
    # define and test coupling
    coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
        section_models)
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Wagner Lifting Line + GXBeam" begin
    # number of elements
    nelem = 2
    # connectivity
    start = 1:nelem
    stop = 2:nelem+1
    # boundary conditions
    displacement = zeros(Bool, 6, nelem+1)
    displacement[:,1] .= true # fixed at one end
    # submodels
    lifting_elements = 1:nelem
    section_models = fill(wagner_liftingline_model(), nelem)
    # define and test coupling
    coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
        section_models)
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

@testset "Peters Lifting Line + GXBeam" begin
    # number of elements
    nelem = 2
    # connectivity
    start = 1:nelem
    stop = 2:nelem+1
    # boundary conditions
    displacement = zeros(Bool, 6, nelem+1)
    displacement[:,1] .= true # fixed at one end
    # submodels
    lifting_elements = 1:nelem
    section_models = fill(peters_liftingline_model(4), nelem)
    # define and test coupling
    coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
        section_models)
    coupling = coupled_model.coupling
    run_coupling_jacobian_tests(coupling)
    run_coupling_io_tests(coupling)
end

# @testset "Steady Lifting Line + Rigid Body" begin
#     coupled_model = liftingline_rigidbody_model(fill(steady_model(), 2))
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# @testset "Quasi-Steady Lifting Line + Rigid Body" begin
#     coupled_model = liftingline_rigidbody_model(fill(quasisteady_model(), 2))
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# @testset "Wagner Lifting Line + Rigid Body" begin
#     coupled_model = liftingline_rigidbody_model(fill(wagner_model(), 2))
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# @testset "Peters Lifting Line + Rigid Body" begin
#     coupled_model = liftingline_rigidbody_model(fill(peters_model(4), 2))
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# @testset "Steady Lifting Line + GXBeam + Rigid Body" begin
#     # number of elements
#     nelem = 2
#     # connectivity
#     start = 1:nelem
#     stop = 2:nelem+1
#     # boundary conditions
#     displacement = zeros(Bool, 6, nelem+1)
#     displacement[:,1] .= true # fixed at one end
#     # submodels
#     lifting_elements = 1:nelem
#     section_models = fill(steady_model(), nelem)
#     # define and test coupling
#     coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
#         section_models)
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# @testset "Quasi-Steady Lifting Line + GXBeam + Rigid Body" begin
#     # number of elements
#     nelem = 2
#     # connectivity
#     start = 1:nelem
#     stop = 2:nelem+1
#     # boundary conditions
#     displacement = zeros(Bool, 6, nelem+1)
#     displacement[:,1] .= true # fixed at one end
#     # submodels
#     lifting_elements = 1:nelem
#     section_models = fill(quasisteady_model(), nelem)
#     # define and test coupling
#     coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
#         section_models)
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# @testset "Wagner Lifting Line + GXBeam + Rigid Body" begin
#     # number of elements
#     nelem = 2
#     # connectivity
#     start = 1:nelem
#     stop = 2:nelem+1
#     # boundary conditions
#     displacement = zeros(Bool, 6, nelem+1)
#     displacement[:,1] .= true # fixed at one end
#     # submodels
#     lifting_elements = 1:nelem
#     section_models = fill(wagner_model(), nelem)
#     # define and test coupling
#     coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
#         section_models)
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# @testset "Peters Lifting Line + GXBeam + Rigid Body" begin
#     # number of elements
#     nelem = 2
#     # connectivity
#     start = 1:nelem
#     stop = 2:nelem+1
#     # boundary conditions
#     displacement = zeros(Bool, 6, nelem+1)
#     displacement[:,1] .= true # fixed at one end
#     # submodels
#     lifting_elements = 1:nelem
#     section_models = fill(peters_model(), nelem)
#     # define and test coupling
#     coupled_model = liftingline_gxbeam_model(start, stop, displacement, lifting_elements, 
#         section_models)
#     coupling = coupled_model.coupling
#     run_coupling_jacobian_tests(coupling)
#     run_coupling_io_tests(coupling)
# end

# # Control Surface Models (3D)

# @testset "LiftingLineFlaps" begin

#     N = 2 # number of lifting line sections
#     NG = 1 # number of control inputs
#     model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     run_model_jacobian_tests(model)
#     run_model_io_tests(model)

# end

# # Controller Models (3D)

# @testset "Trim" begin

#     model = Trim()
#     run_model_jacobian_tests(model)
#     run_model_io_tests(model)

# end

# # --- Three Dimensional Couplings --- #

# @testset "Lifting Line + GXBeam + RigidBody" begin

#     # number of lifting line sections
#     N = 2

#     # point locations
#     x = range(0, 60, length=N+1)
#     y = zero(x)
#     z = zero(x)
#     points = [[x[i],y[i],z[i]] for i = 1:length(x)]

#     # element connectivity
#     start = 1:N
#     stop = 2:N+1

#     # element orientation
#     e1 = [1, 0, 0]
#     e2 = [0, 1, 0]
#     e3 = [0, 0, 1]
#     frame = hcat(e1, e2, e3)

#     # element stiffness properties
#     stiffness = [
#          2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
#          1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
#          6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
#         -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
#         -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
#         -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
#         ]
#     compliance = inv(stiffness)
#     compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

#     # element inertial properties
#     mass = [
#          258.053      0.0        0.0      0.0      7.07839  -71.6871
#            0.0      258.053      0.0     -7.07839  0.0        0.0
#            0.0        0.0      258.053   71.6871   0.0        0.0
#            0.0       -7.07839   71.6871  48.59     0.0        0.0
#            7.07839    0.0        0.0      0.0      2.172      0.0
#          -71.6871     0.0        0.0      0.0      0.0       46.418
#          ]

#     # construct assembly
#     assembly = Assembly(points, start, stop;
#         frames = fill(frame, N),
#         stiffness = fill(stiffness, N),
#         mass = fill(mass, N))

#     # define boundary conditions and applied loads
#     prescribed = Dict(
#             # fixed left side
#             1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
#             # force on right side
#             N+1 => PrescribedConditions(Fz=1e5)
#             )

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{0}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     models = (aerodynamic_model, structural_model, dynamics_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{1}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     models = (aerodynamic_model, structural_model, dynamics_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{2}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     models = (aerodynamic_model, structural_model, dynamics_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Wagner())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     models = (aerodynamic_model, structural_model, dynamics_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Peters{4}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     models = (aerodynamic_model, structural_model, dynamics_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

# end

# @testset "Lifting Line + RigidBody + Lifting Line Flaps" begin

#     # number of lifting line sections
#     N = 2

#     # number of control surface inputs
#     NG = 1

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{0}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{1}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{2}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Wagner())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Peters{4}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

# end

# @testset "Lifting Line + GXBeam + Lifting Line Flaps" begin

#     # number of lifting line sections
#     N = 2

#     # number of control surface inputs
#     NG = 1

#     # point locations
#     x = range(0, 60, length=N+1)
#     y = zero(x)
#     z = zero(x)
#     points = [[x[i],y[i],z[i]] for i = 1:length(x)]

#     # element connectivity
#     start = 1:N
#     stop = 2:N+1

#     # element orientation
#     e1 = [1, 0, 0]
#     e2 = [0, 1, 0]
#     e3 = [0, 0, 1]
#     frame = hcat(e1, e2, e3)

#     # element stiffness properties
#     stiffness = [
#          2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
#          1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
#          6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
#         -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
#         -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
#         -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
#         ]
#     compliance = inv(stiffness)
#     compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

#     # element inertial properties
#     mass = [
#          258.053      0.0        0.0      0.0      7.07839  -71.6871
#            0.0      258.053      0.0     -7.07839  0.0        0.0
#            0.0        0.0      258.053   71.6871   0.0        0.0
#            0.0       -7.07839   71.6871  48.59     0.0        0.0
#            7.07839    0.0        0.0      0.0      2.172      0.0
#          -71.6871     0.0        0.0      0.0      0.0       46.418
#          ]

#     # construct assembly
#     assembly = Assembly(points, start, stop;
#         frames = fill(frame, N),
#         stiffness = fill(stiffness, N),
#         mass = fill(mass, N))

#     # define boundary conditions and applied loads
#     prescribed = Dict(
#             # fixed left side
#             1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
#             # force on right side
#             N+1 => PrescribedConditions(Fz=1e5)
#             )

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{0}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{1}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{2}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Wagner())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Peters{4}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

# end

# @testset "Lifting Line + GXBeam + RigidBody + Lifting Line Flaps" begin

#     # number of lifting line sections
#     N = 2

#     # number of control surface inputs
#     NG = 1

#     # point locations
#     x = range(0, 60, length=N+1)
#     y = zero(x)
#     z = zero(x)
#     points = [[x[i],y[i],z[i]] for i = 1:length(x)]

#     # element connectivity
#     start = 1:N
#     stop = 2:N+1

#     # element orientation
#     e1 = [1, 0, 0]
#     e2 = [0, 1, 0]
#     e3 = [0, 0, 1]
#     frame = hcat(e1, e2, e3)

#     # element stiffness properties
#     stiffness = [
#          2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
#          1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
#          6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
#         -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
#         -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
#         -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
#         ]
#     compliance = inv(stiffness)
#     compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

#     # element inertial properties
#     mass = [
#          258.053      0.0        0.0      0.0      7.07839  -71.6871
#            0.0      258.053      0.0     -7.07839  0.0        0.0
#            0.0        0.0      258.053   71.6871   0.0        0.0
#            0.0       -7.07839   71.6871  48.59     0.0        0.0
#            7.07839    0.0        0.0      0.0      2.172      0.0
#          -71.6871     0.0        0.0      0.0      0.0       46.418
#          ]

#     # construct assembly
#     assembly = Assembly(points, start, stop;
#         frames = fill(frame, N),
#         stiffness = fill(stiffness, N),
#         mass = fill(mass, N))

#     # define boundary conditions and applied loads
#     prescribed = Dict(
#             # fixed left side
#             1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
#             # force on right side
#             N+1 => PrescribedConditions(Fz=1e5)
#             )

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{0}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{1}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{2}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Wagner())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Peters{4}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)
# end

# @testset "Lifting Line + RigidBody + Lifting Line Flaps + Trim" begin

#     # number of lifting line sections
#     N = 2

#     # number of control surface inputs
#     NG = 6

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{0}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{1}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{2}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Wagner())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Peters{4}())
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)
# end

# @testset "Lifting Line + GXBeam + Lifting Line Flaps + Trim" begin

#     # number of lifting line sections
#     N = 2

#     # number of control surface inputs
#     NG = 6

#     # point locations
#     x = range(0, 60, length=N+1)
#     y = zero(x)
#     z = zero(x)
#     points = [[x[i],y[i],z[i]] for i = 1:length(x)]

#     # element connectivity
#     start = 1:N
#     stop = 2:N+1

#     # element orientation
#     e1 = [1, 0, 0]
#     e2 = [0, 1, 0]
#     e3 = [0, 0, 1]
#     frame = hcat(e1, e2, e3)

#     # element stiffness properties
#     stiffness = [
#          2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
#          1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
#          6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
#         -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
#         -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
#         -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
#         ]
#     compliance = inv(stiffness)
#     compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

#     # element inertial properties
#     mass = [
#          258.053      0.0        0.0      0.0      7.07839  -71.6871
#            0.0      258.053      0.0     -7.07839  0.0        0.0
#            0.0        0.0      258.053   71.6871   0.0        0.0
#            0.0       -7.07839   71.6871  48.59     0.0        0.0
#            7.07839    0.0        0.0      0.0      2.172      0.0
#          -71.6871     0.0        0.0      0.0      0.0       46.418
#          ]

#     # construct assembly
#     assembly = Assembly(points, start, stop;
#         frames = fill(frame, N),
#         stiffness = fill(stiffness, N),
#         mass = fill(mass, N))

#     # define boundary conditions and applied loads
#     prescribed = Dict(
#             # fixed left side
#             1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
#             # force on right side
#             N+1 => PrescribedConditions(Fz=1e5)
#             )

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{0}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{1}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{2}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Wagner())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Peters{4}())
#     structural_model = GEBT(assembly, prescribed)
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)
# end

# @testset "Lifting Line + GXBeam + RigidBody + Lifting Line Flaps + Trim" begin

#     # number of lifting line sections
#     N = 2

#     # number of control surface inputs
#     NG = 6

#     # point locations
#     x = range(0, 60, length=N+1)
#     y = zero(x)
#     z = zero(x)
#     points = [[x[i],y[i],z[i]] for i = 1:length(x)]

#     # element connectivity
#     start = 1:N
#     stop = 2:N+1

#     # element orientation
#     e1 = [1, 0, 0]
#     e2 = [0, 1, 0]
#     e3 = [0, 0, 1]
#     frame = hcat(e1, e2, e3)

#     # element stiffness properties
#     stiffness = [
#          2.389e9  1.524e6  6.734e6 -3.382e7 -2.627e7 -4.736e8
#          1.524e6  4.334e8 -3.741e6 -2.935e5  1.527e7  3.835e5
#          6.734e6 -3.741e6  2.743e7 -4.592e5 -6.869e5 -4.742e6
#         -3.382e7 -2.935e5 -4.592e5  2.167e7 -6.279e5  1.430e6
#         -2.627e7  1.527e7 -6.869e5 -6.279e5  1.970e7  1.209e7
#         -4.736e8  3.835e5 -4.742e6  1.430e6  1.209e7  4.406e8
#         ]
#     compliance = inv(stiffness)
#     compliance_entries = [compliance[i,j] for i = 1:6 for j = i:6]

#     # element inertial properties
#     mass = [
#          258.053      0.0        0.0      0.0      7.07839  -71.6871
#            0.0      258.053      0.0     -7.07839  0.0        0.0
#            0.0        0.0      258.053   71.6871   0.0        0.0
#            0.0       -7.07839   71.6871  48.59     0.0        0.0
#            7.07839    0.0        0.0      0.0      2.172      0.0
#          -71.6871     0.0        0.0      0.0      0.0       46.418
#          ]

#     # construct assembly
#     assembly = Assembly(points, start, stop;
#         frames = fill(frame, N),
#         stiffness = fill(stiffness, N),
#         mass = fill(mass, N))

#     # define boundary conditions and applied loads
#     prescribed = Dict(
#             # fixed left side
#             1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
#             # force on right side
#             N+1 => PrescribedConditions(Fz=1e5)
#             )

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{0}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{1}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(QuasiSteady{2}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Wagner())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

#     aerodynamic_model = liftingline_model{N}(Peters{4}())
#     structural_model = GEBT(assembly, prescribed)
#     dynamics_model = RigidBody()
#     surface_model = LiftingLineFlaps{N}(SimpleFlap(), ntuple(i->ones(N), NG))
#     control_model = Trim()
#     models = (aerodynamic_model, structural_model, dynamics_model, surface_model, control_model)
#     run_coupling_jacobian_tests(models...)
#     run_coupling_io_tests(models...)

# end
