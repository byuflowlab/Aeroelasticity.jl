@testset "Wagner" begin

    submodel = Submodel(Wagner())
    
    # test provided jacobians
    submodel_jacobian_tests(submodel)

    # test consistency of input/output functions
    submodel_io_tests(submodel)
end