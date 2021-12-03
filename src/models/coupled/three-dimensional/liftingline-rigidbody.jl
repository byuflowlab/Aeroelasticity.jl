"""
    liftingline_rigidbody_model()

Construct a model by coupling a lifting line aerodynamic model and a rigid body
dynamics model.  This model introduces additional parameters corresponding to
the length, position, and orientation of each lifting line element ``L, p_e,
e_1, e_2, e_3``, followed by the freestream air density ``rho_\\infty``, air speed of sound
``c``, rigid body inertial properties ``m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz``, and the 
additional forces/moments applied on the body ``F_x, F_y, F_z, M_x, M_y, M_z``.

**NOTE: When using this model, the local frame for each lifting line element should be
oriented with the x-axis in the chordwise direction, the y-axis in the spanwise
direction (out the right wing), and the z-axis in the airfoil normal direction**
"""
function liftingline_rigidbody_model(section_models)

    # aerodynamic model
    aero = LiftingLine(section_models)

    # structural model
    dyn = RigidBody()

    # submodels
    submodels = (aero, dyn)

    # construct coupling
    coupling = liftingline_rigidbody_coupling(aero, dyn)

    return CoupledModel(submodels, coupling)
end

# --- Internal Methods --- #

function liftingline_rigidbody_coupling(aero, dyn)

    # extract section models
    section_models = aero.constants.submodels

    # number of aerodynamic sections
    n = length(section_models)

    # number of additional parameters
    npc = 13*n + 15

    # dimensions
    nx = number_of_states(aero) + number_of_states(dyn)
    ny = number_of_inputs(aero) + number_of_inputs(dyn)
    np = number_of_parameters(aero) + number_of_parameters(dyn) + npc

    # coupling function
    g = (dx, x, p, t) -> liftingline_rigidbody_inputs(dx, x, p, t; aero, dyn)

    # jacobians
    ratejac = Nonlinear() # TODO: define rate jacobian
    statejac = Nonlinear() # TODO: define state jacobian
    paramjac = Nonlinear() # TODO: define parameter jacobian
    tgrad = Nonlinear() # TODO: define time gradient
    
    # input/output functions
    setparam = liftingline_rigidbody_setparam!
    sepparam = (p) -> liftingline_rigidbody_sepparam(p; n)

    return Coupling{true}(g, nx, ny, np, npc;
        ratejac = ratejac,
        statejac = statejac,
        paramjac = paramjac,
        tgrad = tgrad,
        setparam = setparam,
        sepparam = sepparam)
end

function liftingline_rigidbody_inputs!(y, dx, x, p, t; aero, dyn)

    # number of aerodynamic sections
    n = length(aero.constants.submodels)

    # get indices of states, inputs, and parameters for each submodel
    submodels = (aero, dyn)
    ixa, ixd = state_indices(submodels)
    iya, iyd = input_indices(submodels)
    ipa, ipd = parameter_indices(submodels)

    # get views of states, inputs, and parameters for each submodel
    dxa, dxd = view(dx, ixa), view(dx, ixd)
    xa, xd = view(x, ixa), view(x, ixd)
    ya, yd = view(y, iya), view(y, iyd)
    pa, pd = view(p, ipa), view(p, ipd) 

    # separate lifting line model rates, states, and parameters
    liftingline_rates = separate_rates(aero, dxa)
    liftingline_states = separate_states(aero, xa)
    liftingline_parameters = separate_parameters(aero, pa)

    # separate rigid body model rates, states, and parameters
    rigidbody_rates = separate_rates(aero, dxd)
    rigidbody_states = separate_states(aero, xd)
    rigidbody_parameters = separate_parameters(aero, pd)

    # separate additional parameters
    np = length(p)
    ipc = np - 13*n - npc : np
    pc = view(p, ipc)
    additional_parameters = liftingline_rigidbody_sepparam(pc; n)

    # initialize loads with body loads
    F = SVector(additional_parameters.Fbx, additional_parameters.Fby, additional_parameters.Fbz)
    M = SVector(additional_parameters.Mbx, additional_parameters.Mby, additional_parameters.Mbz)

    # find environmental conditions
    V, Ω, a, α = rigidbody_environmental_conditions(rigidbody_rates, rigidbody_states)

    # loop through each section
    for i = 1:n

        # get additional parameters for this section
        param = additional_parameters.section_parameters[i]

        # unpack section properties        
        L = param.L
        pe = param.pe
        e1 = param.e1
        e2 = param.e2
        e3 = param.e3

        # unpack section aerodynamic model, rates, states, and parameters
        aerodynamic_model = section_models[i].models[1]
        aerodynamic_rates = liftingline_rates.sections[i]
        aerodynamic_states = liftingline_states.sections[i]
        aerodynamic_parameters = liftingline_parameters.sections[i]

        # determine section velocities/accelerations
        vi, ωi, ai, αi = liftingline_section_velocities(pe, e1, e2, e3, V, Ω, a, α)

        # express section aerodynamic rates, states, and parameters as vectors
        dxai = get_rates(aerodynamic_model; aerodynamic_rates...)
        xai = get_states(aerodynamic_model; aerodynamic_states...)
        pai = get_parameters(aerodynamic_model; aerodynamic_parameters...)

        # express section structural rates, states, and parameters as vectors
        dxsi = SVector(ai..., αi...)
        xsi = SVector(vi..., ωi...)
        psi = SVector(ρ, c)

        # define section rate, state, and parameter vectors
        dxi = SVector(dxai..., dxsi...)
        xi = SVector(xai..., xsi...)
        pi = SVector(pai..., psi...)

        # calculate section inputs
        yi = get_coupling_inputs(section_models[i], dxi, xi, pi, t)

        # separate section inputs
        section_inputs = separate_inputs(section_models[i], yi)

        # separate section aerodynamic and structural inputs 
        aerodynamic_inputs = section_inputs[1]
        structural_inputs = section_inputs[2]

        # save section aerodynamic inputs
        set_inputs!(yai, aerodynamic_model, aerodynamic_inputs)

        # integrate section distributed loads
        Fi = L*structural_inputs.fi
        Mi = L*structural_inputs.mi

        # add section loads to total loads
        F += Fi
        M += cross(pe, Fi) + Mi
    end

    # save rigid body inputs
    set_inputs!(yd, rigidbody_model;
        m = additional_parameters.m, 
        Ixx = additional_parameters.Ixx, 
        Iyy = additional_parameters.Iyy, 
        Izz = additional_parameters.Izz, 
        Ixz = additional_parameters.Ixz, 
        Ixy = additional_parameters.Ixy, 
        Iyz = additional_parameters.Iyz, 
        Fx = F[1], 
        Fy = F[2], 
        Fz = F[3], 
        Mx = M[1], 
        My = M[2], 
        Mz = M[3])

    return y
end

function liftingline_rigidbody_setparam!(p; section_parameters, rho, c,
    m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fbx, Fby, Fbz, Mbx, Mby, Mbz)

    n = length(section_parameters)

    for i = 1:n
        offset = 13*(i-1)
        p[offset+1] = section_parameters[i].L
        p[offset+2:offset+4] = section_parameters[i].pe
        p[offset+5:offset+7] = section_parameters[i].e1
        p[offset+8:offset+10] = section_parameters[i].e2
        p[offset+11:offset+13] = section_parameters[i].e3
    end

    offset = 13*(n-1)
    p[offset+1] = rho
    p[offset+2] = c
    p[offset+3] = m
    p[offset+4] = Ixx
    p[offset+5] = Iyy
    p[offset+6] = Izz
    p[offset+7] = Ixz
    p[offset+8] = Ixy
    p[offset+9] = Iyz
    p[offset+10] = Fbx
    p[offset+11] = Fby
    p[offset+12] = Fbz
    p[offset+13] = Mbx
    p[offset+14] = Mby
    p[offset+15] = Mbz

    return p
end

liftingline_rigidbody_sepparam(p; n) = (
    section_parameters = ntuple(i -> (
        L = p[13*(i-1)+1], # element length
        pe = view(p, 13*(i-1)+2 : 13*(i-1)+4), # local frame position
        e1 = view(p, 13*(i-1)+5 : 13*(i-1)+7), # local frame x-axis
        e2 = view(p, 13*(i-1)+8 : 13*(i-1)+10), # local frame y-axis
        e3 = view(p, 13*(i-1)+11 : 13*(i-1)+13), # local frame z-axis
        ), n),
    rho = p[13*n+1],
    c = p[13*n+2],
    m = p[13*n+3],
    Ixx = p[13*n+4],
    Iyy = p[13*n+5],
    Izz = p[13*n+6],
    Ixz = p[13*n+7],
    Ixy = p[13*n+8],
    Iyz = p[13*n+9],
    Fbx = p[13*n+10],
    Fby = p[13*n+11],
    Fbz = p[13*n+12],
    Mbx = p[13*n+13],
    Mby = p[13*n+14],
    Mbz = p[13*n+15],
    )


