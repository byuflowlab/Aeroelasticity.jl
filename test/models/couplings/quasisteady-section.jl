@testset "QuasiSteady + Section" begin
    
    # define models
    models = (QuasiSteady(), Section())

    # create coupling
    coupling = Coupling(models)

    # test provided jacobians
    coupling_jacobian_tests(coupling)
    
    # test consistency of input/output functions
    coupling_io_tests(coupling)
end