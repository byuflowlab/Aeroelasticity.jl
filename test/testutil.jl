function submodel_jacobian_tests(submodel;
    dx = rand(Aeroelasticity.number_of_states(submodel)),
    x = rand(Aeroelasticity.number_of_states(submodel)),
    y = rand(Aeroelasticity.number_of_inputs(submodel)),
    p = rand(Aeroelasticity.number_of_parameters(submodel)),
    t = rand(),
    atol = sqrt(eps()),
    norm = (x) -> LinearAlgebra.norm(x, Inf))

    # rate jacobian test
    f = (dx) -> Aeroelasticity.get_residual(submodel, dx, x, y, p, t)
    J = Aeroelasticity.get_rate_jacobian(submodel, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, dx)
    @test isapprox(J, J_fd; atol, norm)

    # state jacobian test
    f = (x) -> Aeroelasticity.get_residual(submodel, dx, x, y, p, t)
    J = Aeroelasticity.get_state_jacobian(submodel, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, x)
    @test isapprox(J, J_fd; atol, norm)

    # input jacobian test
    f = (y) -> Aeroelasticity.get_residual(submodel, dx, x, y, p, t)
    J = Aeroelasticity.get_input_jacobian(submodel, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, y)
    @test isapprox(J, J_fd; atol, norm)

    # parameter jacobian test
    f = (p) -> Aeroelasticity.get_residual(submodel, dx, x, y, p, t)
    J = Aeroelasticity.get_parameter_jacobian(submodel, dx, x, y, p, t)
    J_fd = ForwardDiff.jacobian(f, p)
    @test isapprox(J, J_fd; atol, norm)

    # time derivative test
    f = (t) -> Aeroelasticity.get_residual(submodel, dx, x, y, p, t)
    dT = Aeroelasticity.get_time_gradient(submodel, dx, x, y, p, t)
    dT_fd = ForwardDiff.derivative(f, t)
    @test isapprox(dT, dT_fd; atol, norm)

    return nothing
end

function coupling_jacobian_tests(coupling;
    dx = rand(Aeroelasticity.number_of_states(coupling)),
    x = rand(Aeroelasticity.number_of_states(coupling)),
    p = rand(Aeroelasticity.number_of_parameters(coupling)),
    t = rand(),
    atol = sqrt(eps()),
    norm = (x)->norm(x, Inf))

    # rate jacobian test
    f = (dx) -> Aeroelasticity.get_coupling_inputs(coupling, dx, x, p, t)
    J = Aeroelasticity.get_coupling_rate_jacobian(coupling, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, dx)
    @test isapprox(J, J_fd; atol, norm)

    # state jacobian test
    f = (x) -> Aeroelasticity.get_coupling_inputs(coupling, dx, x, p, t)
    J = Aeroelasticity.get_coupling_state_jacobian(coupling, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, x)
    @test isapprox(J, J_fd; atol, norm)

    # parameter jacobian test
    f = (p) -> Aeroelasticity.get_coupling_inputs(coupling, dx, x, p, t)
    J = Aeroelasticity.get_coupling_parameter_jacobian(coupling, dx, x, p, t)
    J_fd = ForwardDiff.jacobian(f, p)
    @test isapprox(J, J_fd; atol, norm)

    # time derivative test
    f = (t) -> Aeroelasticity.get_coupling_inputs(coupling, dx, x, p, t)
    dT = Aeroelasticity.get_coupling_time_gradient(coupling, dx, x, p, t)
    dT_fd = ForwardDiff.derivative(f, t)
    @test isapprox(dT, dT_fd; atol, norm)

    return nothing
end

function submodel_io_tests(model; atol = sqrt(eps()), norm = (x)->norm(x, Inf))

    # test separating and recombining state variables
    x = rand(Aeroelasticity.number_of_states(model))
    states = Aeroelasticity.separate_states(model, x)
    new_x = Aeroelasticity.set_states(model; states...)
    @test isapprox(new_x, x; atol, norm)

    # test separating and recombining inputs
    y = rand(Aeroelasticity.number_of_inputs(model))
    inputs = Aeroelasticity.separate_inputs(model, y)
    new_y = Aeroelasticity.set_inputs(model; inputs...)
    @test isapprox(new_y, y; atol, norm)

    # test separating and recombining parameters
    p = rand(Aeroelasticity.number_of_parameters(model))
    parameters = Aeroelasticity.separate_parameters(model, p)
    new_p = Aeroelasticity.set_parameters(model; parameters...)
    @test isapprox(new_p, p; atol, norm)

    return nothing
end

function coupling_io_tests(coupling; atol = sqrt(eps()), norm = (x)->norm(x, Inf))

    # test separating and recombining additional parameters
    pc = rand(Aeroelasticity.number_of_additional_parameters(coupling))
    parameters = Aeroelasticity.separate_parameters(coupling, pc)
    new_pc = Aeroelasticity.set_parameters(coupling; parameters...)
    @test isapprox(new_pc, pc; atol, norm)

    return nothing
end