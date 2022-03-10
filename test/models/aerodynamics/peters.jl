@testset "Peters' Finite State" begin

    submodel = Submodel(Peters{4}())

    # test provided jacobians
    submodel_jacobian_tests(submodel)

    # test consistency of input/output functions
    submodel_io_tests(submodel)
end
