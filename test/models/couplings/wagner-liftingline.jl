@testset "Wagner + LiftingLineSection" begin
    
    # define models
    models = (Wagner(), LiftingLineSection())

    # create coupling
    coupling = Coupling(models)

    # test provided jacobians
    coupling_jacobian_tests(coupling)

    # test consistency of input/output functions
    coupling_io_tests(coupling)
end