@testset "LiftingLineSection" begin

    submodel = Submodel(LiftingLineSection())
    
    # test consistency of input/output functions
    submodel_io_tests(submodel)
end