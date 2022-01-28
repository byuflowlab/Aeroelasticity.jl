@testset "Lifting Line" begin

    section_models = ntuple(Steady(), 2)
    submodel = Submodel(LiftingLine(section_models))
    model_jacobian_tests(submodel)
    model_io_tests(submodel)

    section_models = ntuple(QuasiSteady(), 2)
    submodel = Submodel(LiftingLine(section_models))
    model_jacobian_tests(submodel)
    model_io_tests(submodel)

    section_models = ntuple(Wagner(), 2)
    submodel = Submodel(LiftingLine(section_models))
    model_jacobian_tests(submodel)
    model_io_tests(submodel)

    section_models = ntuple(Peters(), 2)
    submodel = Submodel(LiftingLine(section_models))
    model_jacobian_tests(submodel)
    model_io_tests(submodel)

end