@testset "LiftingLine + GXBeamAssembly" begin

    # number of elements
    nelem = 2

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
    prescribed_conditions = Dict(
            # fixed left side
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0),
            # force on right side
            nelem+1 => PrescribedConditions(Fz=1e5)
            )

    section_models = fill(Steady(), nelem)
    models = (LiftingLine(section_models), GXBeamAssembly(assembly, prescribed_conditions))
    coupling = Coupling(models)
    coupling_jacobian_tests(coupling)
    coupling_io_tests(coupling)

    section_models = fill(QuasiSteady(), nelem)
    models = (LiftingLine(section_models), GXBeamAssembly(assembly, prescribed_conditions))
    coupling = Coupling(models)
    coupling_jacobian_tests(coupling)
    coupling_io_tests(coupling)

    section_models = fill(Wagner(), nelem)
    models = (LiftingLine(section_models), GXBeamAssembly(assembly, prescribed_conditions))
    coupling = Coupling(models)
    coupling_jacobian_tests(coupling)
    coupling_io_tests(coupling)

    section_models = fill(Peters{4}(), nelem)
    models = (LiftingLine(section_models), GXBeamAssembly(assembly, prescribed_conditions))
    coupling = Coupling(models)
    coupling_jacobian_tests(coupling)
    coupling_io_tests(coupling)
        
end