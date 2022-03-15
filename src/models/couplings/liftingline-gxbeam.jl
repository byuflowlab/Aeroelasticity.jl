# --- Coupling Model Creation --- #

"""
    Coupling(::LiftingLine, ::GXBeamAssembly; kwargs...)

Coupling model for coupling a lifting line aerodynamic model and a geometrically exact beam 
element assembly model.  This model introduces additional parameters corresponding to
the external point loads ``F_{x,i}, F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i}, M_{z,i}`` or 
displacements ``u_{x,i}, u_{y,i}, u_{z,i}, \\theta_{x,i}, \\theta_{y,i}, \\theta_{z,i}`` 
applied to each node, followed by the distributed loads ``f_{x,i}, f_{y,i}, f_{z,i}, 
m_{x,i}, m_{y,i}, m_{z,i}`` applied to each beam element (excluding aerodynamic loads), 
followed by the properties of point masses attached to each beam element ``m, p, I_{11}, 
I_{22}, I_{33}, I_{12}, I_{13}, I_{23}``, followed by the gravity vector, linear velocity, 
angular velocity, linear acceleration, and angular acceleration of the system, followed by 
the freestream air density ``\\rho_\\infty`` and air speed of sound ``c``.

# Keyword Arguments
 - `lifting_elements`: Beam element index corresponding to each aerodynamic section.
    By default, aerodynamic models are assigned to beam elements in order

**NOTE: This model assumes that each beam element is oriented with the x-axis along the 
beam's axis, the y-axis forward (into the freestream), and the z-axis normal**
"""
function Coupling(models::Tuple{LiftingLine, GXBeamAssembly}, submodels; 
    lifting_elements = 1:length(models[1].section_models))

    liftingline, gxbeam = models

    npoint = length(gxbeam.icol_point)
    nelem = length(gxbeam.icol_elem)
    displacement = gxbeam.displacement

    # construct a coupled model for each section
    section_models = broadcast(liftingline.section_models) do aerodynamic_model
        assemble_model(; 
            aerodynamic_model = aerodynamic_model,
            structural_model = LiftingLineSection())
    end
    
    # construct the coupling function
    g = (y, dx, x, p, t) -> liftingline_gxbeam_inputs!(y, dx, x, p, t; models, submodels,
        section_models, lifting_elements)
    
    # define the number of additional parameters
    npc = 6*npoint + 6*nelem + 10*nelem + 17

    # define the number of states, inputs, and parameters
    nx = sum(number_of_states.(submodels))
    ny = sum(number_of_inputs.(submodels))
    np = sum(number_of_parameters.(submodels)) + npc

    # define the jacobians
    ratejac = Nonlinear() # TODO: define rate jacobian
    statejac = Nonlinear() # TODO: define state jacobian
    paramjac = Nonlinear() # TODO: define parameter jacobian
    tgrad = Nonlinear() # TODO: define time gradient
    
    # define the input/output functions
    setparam = (p; kwargs...) -> liftingline_gxbeam_setparam!(p, npoint, nelem, displacement; kwargs...)
    sepparam = (p) -> liftingline_gxbeam_sepparam(p, npoint, nelem, displacement)

    return Coupling{true}(g, nx, ny, np, npc;
        ratejac = ratejac,
        statejac = statejac,
        paramjac = paramjac,
        tgrad = tgrad,
        setparam = setparam,
        sepparam = sepparam)
end

# coupling function
function liftingline_gxbeam_inputs!(y, dx, x, p, t; models, submodels, section_models,
    lifting_elements)

    liftingline, gxbeam = models
    aero, stru = submodels

    # get indices of states, inputs, and parameters for each submodel
    ixa, ixs = state_indices(submodels)
    iya, iys = input_indices(submodels)
    ipa, ips = parameter_indices(submodels)
    ipc = sum(number_of_parameters.(submodels)) + 1 : length(p)

    # get views of states, inputs, and parameters for each submodel
    dxa, dxs = view(dx, ixa), view(dx, ixs)
    xa, xs = view(x, ixa), view(x, ixs)
    ya, ys = view(y, iya), view(y, iys)
    pa, ps, pc = view(p, ipa), view(p, ips), view(p, ipc) 

    # get indices of states, inputs, and parameters for each lifting line section model
    section_aerodynamic_submodels = getindex.(getproperty.(section_models, :submodels), 1)
    ixas = state_indices(section_aerodynamic_submodels)
    iyas = input_indices(section_aerodynamic_submodels)
    ipas = parameter_indices(section_aerodynamic_submodels)

    # get views of states, inputs, and parameters for each lifting line section model
    dxas = view.(Ref(dxa), ixas)
    xas = view.(Ref(xa), ixas)
    yas = view.(Ref(ya), iyas)
    pas = view.(Ref(pa), ipas) 

    # separate gxbeam model rates, states, and parameters
    gxbeam_rates = separate_rates(stru, dxs)
    gxbeam_states = separate_states(stru, xs)
    gxbeam_parameters = separate_parameters(stru, ps)

    # separate additional parameters
    nelem = length(gxbeam_parameters.assembly.elements)
    npoint = length(gxbeam_parameters.assembly.points)
    additional_parameters = liftingline_gxbeam_sepparam(pc, npoint, nelem, gxbeam.displacement)

    # unpack element rate variables
    @unpack u_e, theta_e, F_e, M_e, V_e, Omega_e = gxbeam_rates 
    du_e, dtheta_e, dF_e, dM_e, dV_e, dOmega_e = u_e, theta_e, F_e, M_e, V_e, Omega_e

    # unpack element state variables
    @unpack u_e, theta_e, F_e, M_e, V_e, Omega_e = gxbeam_states 

    # unpack addition parameters
    @unpack a, alpha, rho, c, f_d, m_d = additional_parameters
    
    # initialize distributed loads
    f_d = SVector{3,eltype(y)}.(f_d)
    m_d = SVector{3,eltype(y)}.(m_d)

    # calculate aerodynamic inputs and loads for each lifting element
    for (i, ie) in enumerate(lifting_elements)

        # extract element
        elem = gxbeam_parameters.assembly.elements[ie]

        # get element velocity and acceleration (in the element frame)
        vi, ωi = gxbeam_element_velocities(V_e[ie], Omega_e[ie])
        ai, αi = gxbeam_element_accelerations(elem, u_e[ie], theta_e[ie], dV_e[ie], 
            dOmega_e[ie], a, alpha)

        # transform velocities and accelerations to the aerodynamic frame
        local_to_aero = @SMatrix [0 -1 0; 1 0 0; 0 0 1]
        vi, ωi = -local_to_aero*vi, local_to_aero*ωi
        ai, αi = -local_to_aero*ai, local_to_aero*αi

        # calculate lifting line section inputs
        yi = liftingline_section_inputs(section_models[i], dxas[i], xas[i], pas[i], vi, ωi, 
            ai, αi, rho, c, t)

        # separate lifting line section inputs
        aerodynamic_inputs, structural_inputs = separate_inputs(section_models[i], yi)

        # save aerodynamic inputs
        set_inputs!(yas[i], section_models[i].submodels[1]; aerodynamic_inputs...)

        # transform loads to the body frame
        body_to_local = gxbeam_body_to_element(elem, theta_e[ie])
        aero_to_body = (body_to_local*local_to_aero)'

        # add distributed loads from aerodynamics to prescribed distributed loads
        f_d[ie] += aero_to_body*structural_inputs.f
        m_d[ie] += aero_to_body*structural_inputs.m
    end

    # save structural inputs
    set_inputs!(ys, stru; 
        u_p = additional_parameters.u_p, 
        theta_p = additional_parameters.theta_p, 
        F_p = additional_parameters.F_p, 
        M_p = additional_parameters.M_p,
        f_d = f_d,
        m_d = m_d,
        m_m = additional_parameters.m_m,
        p_m = additional_parameters.p_m,
        I_m = additional_parameters.I_m,
        gvec = additional_parameters.gvec,
        v = additional_parameters.v,
        omega = additional_parameters.omega,
        a = additional_parameters.a,
        alpha = additional_parameters.alpha)
        
    return y
end

# convenience function for defining this model's input vector
function liftingline_gxbeam_setparam!(p, np, ne, displacement; 
    u_p = nothing, theta_p = nothing, F_p = nothing, M_p = nothing, 
    f_d = nothing, m_d = nothing, 
    m_m = nothing, p_m = nothing, I_m = nothing,
    gvec = nothing, v = nothing, omega = nothing, a = nothing, alpha = nothing,
    rho = nothing, c = nothing)

    for ip = 1:np
        prescribed_forces = SVector{6}(view(displacement, :, ip)) .== false
        for i = 1:3
            if prescribed_forces[i]
                if !isnothing(F_p)
                    p[6*(ip-1) + i] = F_p[ip][i]
                end
            else
                if !isnothing(u_p)
                    p[6*(ip-1) + i] = u_p[ip][i]
                end
            end
        end
        for i = 1:3
            if prescribed_forces[i]
                if !isnothing(M_p)
                    p[6*(ip-1) + 3 + i] = M_p[ip][i]
                end
            else
                if !isnothing(theta_p)
                    p[6*(ip-1) + 3 + i] = theta_p[ip][i]
                end
            end
        end
    end

    for ie = 1:ne
        if !isnothing(f_d)
            for i = 1:3
                p[6*np + 6*(ie-1) + i] = f_d[ie][i]
            end
        end
        if !isnothing(m_d)
            for i = 1:3
                p[6*np + 6*(ie-1) + 3 + i] = m_d[ie][i]
            end
        end
    end

    for ie = 1:ne
        if !isnothing(m_m)
            p[6*np + 6*ne + 10*(ie-1) + 1] = m_m[ie]
        end
        if !isnothing(p_m)
            p[6*np + 6*ne + 10*(ie-1) + 2] = p_m[ie][1]
            p[6*np + 6*ne + 10*(ie-1) + 3] = p_m[ie][2]
            p[6*np + 6*ne + 10*(ie-1) + 4] = p_m[ie][3]
        end
        if !isnothing(I_m)
            p[6*np + 6*ne + 10*(ie-1) + 5] = I_m[ie][1,1]
            p[6*np + 6*ne + 10*(ie-1) + 6] = I_m[ie][2,2]
            p[6*np + 6*ne + 10*(ie-1) + 7] = I_m[ie][3,3]
            p[6*np + 6*ne + 10*(ie-1) + 8] = I_m[ie][1,2]
            p[6*np + 6*ne + 10*(ie-1) + 9] = I_m[ie][1,3]
            p[6*np + 6*ne + 10*(ie-1) + 10] = I_m[ie][2,3]
        end
    end

    if !isnothing(gvec)
        p[6*np + 6*ne + 10*ne + 1] = gvec[1]
        p[6*np + 6*ne + 10*ne + 2] = gvec[2]
        p[6*np + 6*ne + 10*ne + 3] = gvec[3]
    end

    if !isnothing(v)
        p[6*np + 6*ne + 10*ne + 4] = v[1]
        p[6*np + 6*ne + 10*ne + 5] = v[2]
        p[6*np + 6*ne + 10*ne + 6] = v[3]
    end

    if !isnothing(omega)
        p[6*np + 6*ne + 10*ne + 7] = omega[1]
        p[6*np + 6*ne + 10*ne + 8] = omega[2]
        p[6*np + 6*ne + 10*ne + 9] = omega[3]
    end

    if !isnothing(a)
        p[6*np + 6*ne + 10*ne + 10] = a[1]
        p[6*np + 6*ne + 10*ne + 11] = a[2]
        p[6*np + 6*ne + 10*ne + 12] = a[3]
    end

    if !isnothing(alpha)
        p[6*np + 6*ne + 10*ne + 13] = alpha[1]
        p[6*np + 6*ne + 10*ne + 14] = alpha[2]
        p[6*np + 6*ne + 10*ne + 15] = alpha[3]
    end

    if !isnothing(rho)
        p[6*np + 6*ne + 10*ne + 16] = rho
    end

    if !isnothing(c)
        p[6*np + 6*ne + 10*ne + 17] = c
    end

    return p
end

# convenience function for separating this model's parameter vector
function liftingline_gxbeam_sepparam(p, np, ne, displacement)

    TF = eltype(p)

    # separate prescribed conditions
    u_p = Vector{SVector{3,TF}}(undef, np)
    theta_p = Vector{SVector{3,TF}}(undef, np)
    F_p = Vector{SVector{3,TF}}(undef, np)
    M_p = Vector{SVector{3,TF}}(undef, np)
    for ip = 1:np
        prescribed_forces = SVector{6}(view(displacement, :, ip)) .== false
        F_p[ip] = SVector(
            ifelse(prescribed_forces[1], p[6*(ip-1) + 1], NaN),
            ifelse(prescribed_forces[2], p[6*(ip-1) + 2], NaN),
            ifelse(prescribed_forces[3], p[6*(ip-1) + 3], NaN))
        M_p[ip] = SVector(
            ifelse(prescribed_forces[4], p[6*(ip-1) + 4], NaN),
            ifelse(prescribed_forces[5], p[6*(ip-1) + 5], NaN),
            ifelse(prescribed_forces[6], p[6*(ip-1) + 6], NaN))
        u_p[ip] = SVector(
            ifelse(prescribed_forces[1], NaN, p[6*(ip-1) + 1]),
            ifelse(prescribed_forces[2], NaN, p[6*(ip-1) + 2]),
            ifelse(prescribed_forces[3], NaN, p[6*(ip-1) + 3]))
        theta_p[ip] = SVector(
            ifelse(prescribed_forces[4], NaN, p[6*(ip-1) + 4]),
            ifelse(prescribed_forces[5], NaN, p[6*(ip-1) + 5]),
            ifelse(prescribed_forces[6], NaN, p[6*(ip-1) + 6]))
    end

    # separate distributed loads
    f_d = Vector{SVector{3,TF}}(undef, ne)
    m_d = Vector{SVector{3,TF}}(undef, ne)
    for ie = 1:ne
        f_d[ie] = SVector(
            p[6*np + 6*(ie-1) + 1], 
            p[6*np + 6*(ie-1) + 2], 
            p[6*np + 6*(ie-1) + 3])
        m_d[ie] = SVector(
            p[6*np + 6*(ie-1) + 4], 
            p[6*np + 6*(ie-1) + 5], 
            p[6*np + 6*(ie-1) + 6])
    end

    # separate point masses
    m_m = Vector{TF}(undef, ne)
    p_m = Vector{SVector{3,TF}}(undef, ne)
    I_m = Vector{SMatrix{3,3,TF}}(undef, ne)
    for ie = 1:ne
        m_m[ie] = p[6*np + 6*ne + 10*(ie-1) + 1]
        p_m[ie] = SVector(
            p[6*np + 6*ne + 10*(ie-1) + 2], 
            p[6*np + 6*ne + 10*(ie-1) + 3], 
            p[6*np + 6*ne + 10*(ie-1) + 4])
        I11 = p[6*np + 6*ne + 10*(ie-1) + 5]
        I22 = p[6*np + 6*ne + 10*(ie-1) + 6]
        I33 = p[6*np + 6*ne + 10*(ie-1) + 7]
        I12 = p[6*np + 6*ne + 10*(ie-1) + 8]
        I13 = p[6*np + 6*ne + 10*(ie-1) + 9]
        I23 = p[6*np + 6*ne + 10*(ie-1) + 10]
        I_m[ie] = @SMatrix [I11 I12 I13; I12 I22 I23; I13 I23 I33]
    end

    # gravity
    gvec =  view(p, 6*np + 6*ne + 10*ne +  1 : 6*np + 6*ne + 10*ne +  3)
    
    # linear/angular velocity
    v =     view(p, 6*np + 6*ne + 10*ne +  4 : 6*np + 6*ne + 10*ne +  6)
    omega = view(p, 6*np + 6*ne + 10*ne +  7 : 6*np + 6*ne + 10*ne +  9)
    
    # linear/angular acceleration
    a =     view(p, 6*np + 6*ne + 10*ne + 10 : 6*np + 6*ne + 10*ne + 12)
    alpha = view(p, 6*np + 6*ne + 10*ne + 13 : 6*np + 6*ne + 10*ne + 15)

    # freestream air density
    rho = p[6*np + 6*ne + 10*ne + 16]

    # air speed of sound
    c = p[6*np + 6*ne + 10*ne + 17]

    return (u_p=u_p, theta_p=theta_p, F_p=F_p, M_p=M_p, f_d=f_d, m_d=m_d, m_m=m_m, p_m=p_m, 
        I_m=I_m, gvec=gvec, v=v, omega=omega, a=a, alpha=alpha, rho=rho, c=c)
end
