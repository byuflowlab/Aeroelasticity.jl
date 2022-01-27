@testset "Peters' Finite State + Section" begin

    coupling = Coupling(Peters{4}(), Section())

    # test provided jacobians
    coupling_jacobian_tests(coupling)

    # test consistency of input/output functions
    coupling_io_tests(coupling)
end