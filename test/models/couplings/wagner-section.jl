@testset "Wagner + Section" begin
    
    coupling = Coupling(Wagner(), Section())
    
    # test provided jacobians
    coupling_jacobian_tests(coupling)

    # test consistency of input/output functions
    coupling_io_tests(coupling)
end