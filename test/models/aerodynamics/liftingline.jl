@testset "Lifting Line" begin

    section_models = fill(Steady(), 2)
    submodel = Submodel(LiftingLine(section_models))
    submodel_jacobian_tests(submodel)
    submodel_io_tests(submodel)

    section_models = fill(QuasiSteady(), 2)
    submodel = Submodel(LiftingLine(section_models))
    submodel_jacobian_tests(submodel)
    submodel_io_tests(submodel)

    section_models = fill(Wagner(), 2)
    submodel = Submodel(LiftingLine(section_models))
    submodel_jacobian_tests(submodel)
    submodel_io_tests(submodel)

    section_models = fill(Peters{4}(), 2)
    submodel = Submodel(LiftingLine(section_models))
    submodel_jacobian_tests(submodel)
    submodel_io_tests(submodel)

end