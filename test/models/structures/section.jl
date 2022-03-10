@testset "Typical Section" begin

    submodel = Submodel(Section())
    
    # test provided jacobians
    submodel_jacobian_tests(submodel)
    
    # test consistency of input/output functions
    submodel_io_tests(submodel)
end