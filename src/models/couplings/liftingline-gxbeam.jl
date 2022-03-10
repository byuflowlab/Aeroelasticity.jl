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
function Coupling(aero::LiftingLine, stru::GXBeamAssembly; 
    lifting_elements = 1:length(aero.submodels))

    # section models
    section_models = aero.section_models

    # number of points and elements
    npoint = length(gxbeam.constants.icol_point)
    nelem = length(gxbeam.constants.icol_elem)

    # displacement constraints
    displacement = gxbeam.constants.displacement

    # number of additional parameters
    npc = 6*npoint + 6*nelem + 10*nelem + 17

    # dimensions
    nx = number_of_states(liftingline) + number_of_states(gxbeam)
    ny = number_of_inputs(liftingline) + number_of_inputs(gxbeam)
    np = number_of_parameters(liftingline) + number_of_parameters(gxbeam) + npc

    # coupling function
    g = (y, dx, x, p, t) -> liftingline_gxbeam_inputs!(y, dx, x, p, t; 
        aero, stru, lifting_elements, section_models)

    # jacobians
    ratejac = Nonlinear() # TODO: define rate jacobian
    statejac = Nonlinear() # TODO: define state jacobian
    paramjac = Nonlinear() # TODO: define parameter jacobian
    tgrad = Nonlinear() # TODO: define time gradient
    
    # input/output functions
    setparam = (p; kwargs...) -> liftingline_gxbeam_setparam!(p, displacement, 
        npoint, nelem; kwargs...)
    sepparam = (p) -> liftingline_gxbeam_sepparam(p, displacement, npoint, nelem)

    return Coupling{true}(g, nx, ny, np, npc;
        ratejac = ratejac,
        statejac = statejac,
        paramjac = paramjac,
        tgrad = tgrad,
        setparam = setparam,
        sepparam = sepparam)
end

function liftingline_gxbeam_inputs!(y, dx, x, p, t; 
    liftingline_submodel, gxbeam_submodel, lifting_elements, section_models)

    # get indices of states, inputs, and parameters for each submodel
    submodels = (liftingline, gxbeam)
    ixa, ixs = state_indices(submodels)
    iya, iys = input_indices(submodels)
    ipa, ips = parameter_indices(submodels)

    # get views of states, inputs, and parameters for each submodel
    dxa, dxs = view(dx, ixa), view(dx, ixs)
    xa, xs = view(x, ixa), view(x, ixs)
    ya, ys = view(y, iya), view(y, iys)
    pa, ps = view(p, ipa), view(p, ips) 

    # get indices of states, inputs, and parameters for each lifting line section model
    ixas = state_indices(section_models)
    iyas = input_indices(section_models)
    ipas = parameter_indices(section_models)

    # get views of states, inputs, and parameters for each lifting line section model
    dxas = view.(Ref(dxa), ixas)
    xas = view.(Ref(xa), ixas)
    yas = view.(Ref(ya), iyas)
    pas = view.(Ref(pa), ipas) 

    # separate gxbeam model rates, states, and parameters
    gxbeam_rates = separate_rates(gxbeam, dxs)
    gxbeam_states = separate_states(gxbeam, xs)
    gxbeam_parameters = separate_parameters(gxbeam, ps)

    # separate additional parameters
    np = length(p)
    npc = np - sum(number_of_parameters.(submodels))
    ipc = np - npc : np
    pc = view(p, ipc)
    additional_parameters = liftingline_gxbeam_sepparam(pc, gxbeam.constants.displacement, 
        length(gxbeam.constants.icol_point), length(gxbeam.constants.icol_elem))

    # initialize distributed loads
    f_d = SVector{3,eltype(y)}.(additional_parameters.f_d)
    m_d = SVector{3,eltype(y)}.(additional_parameters.m_d)

    # calculate aerodynamic inputs and loads
    for (i, ielem) in enumerate(lifting_elements)

        # unpack element properties
        elem = gxbeam_parameters.assembly.elements[ielem]

        # transformation from body to element frame
        C1 = gxbeam_body_to_element(elem, gxbeam_states.theta_e[ielem])

        # element velocities and accelerations (in the element frame)
        vi, ωi = gxbeam_element_velocities(gxbeam_states.V_e[ielem], gxbeam_states.Omega_e[ielem])
        ai, αi = gxbeam_element_accelerations(elem, gxbeam_states.u_e[ielem], 
            gxbeam_states.theta_e[ielem], gxbeam_rates.V_e[ielem], gxbeam_rates.Omega_e[ielem], 
            additional_parameters.a, additional_parameters.alpha)

        # transformation from element to local aerodynamic frame
        C2 = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # freestream velocities and accelerations (in the local aerodynamic frame)
        vi, ωi = -C2*vi, C2*ωi
        ai, αi = -C2*ai, C2*αi

        # calculate lifting line inputs
        yi = liftingline_section_inputs(section_models[i], dxas[i], xas[i], pas[i], 
            vi, ωi, ai, αi, additional_parameters.rho, additional_parameters.c, t)

        # separate section aerodynamic and structural inputs
        aerodynamic_inputs, structural_inputs = separate_inputs(section_models[i], yi)

        # save section aerodynamic inputs
        set_inputs!(yas[i], liftingline_section_models[i]; aerodynamic_inputs...)

        # transformation from local aerodynamic frame to body frame
        C3 = C1'*C2' 

        # add aerodynamic loads to prescribed distributed loads
        f_d[ielem] += C3*structural_inputs.f
        m_d[ielem] += C3*structural_inputs.m
    end

    # save gxbeam inputs
    set_inputs!(ys, gxbeam; 
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
function liftingline_gxbeam_setparam!(p, displacement, np, ne; 
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

# convenience function for separating this model's input vector
function liftingline_gxbeam_sepparam(p, displacement, np, ne)

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
