@testset "QuasiSteady + Section" begin
    
    coupling = Coupling(QuasiSteady(), Section())
    
    # test provided jacobians
    coupling_jacobian_tests(coupling)
    
    # test consistency of input/output functions
    coupling_io_tests(coupling)
end