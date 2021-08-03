"""
    couple_models(aero::LiftingLine, stru::GEBT)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a geometrically exact beam theory model.  This model introduces additional
parameters corresponding to the freestream velocity components ``V_x, V_y, V_z``,
followed by the air density ``\\rho``, followed by the external loads ``F_{x,i},
F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i}, M_{z,i}`` or displacements ``u_{x,i},
u_{y,i}, u_{z,i}, \\theta_{x,i}, \\theta_{y,i}, \\theta_{z,i}`` for each node,
followed by the constant distributed loads ``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i},
m_{y,i}, m_{z,i}`` for each beam element (excluding aerodynamic loads), followed
by the body frame linear and angular velocity ``u, v, w, p, q, r``.

**NOTE: When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface**
"""
couple_models(aero::LiftingLine, stru::GEBT) = (aero, stru)

# --- traits --- #

function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function number_of_parameters(aero::LiftingLine, stru::GEBT)
    return 4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 6
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npadd = number_of_parameters(aero, stru) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    ipadd = npa + nps + 1 : npa + nps + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    us = view(u, ius) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # global parameters
    Vx = padd[1]
    Vy = padd[2]
    Vz = padd[3]
    ρ = padd[4]
    ur = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 1]
    vr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 2]
    wr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 3]
    pr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 4]
    qr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 5]
    rr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 6]

    # body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - V

    # construct assembly from structural parameters
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # save prescribed point loads/displacements
    for ip = 1:npoint
        yoff = nya + 6*(ip-1)
        poff = npa + nps + 4 + 6*(ip-1)
        y[yoff+1] = p[poff+1]
        y[yoff+2] = p[poff+2]
        y[yoff+3] = p[poff+3]
        y[yoff+4] = p[poff+4]
        y[yoff+5] = p[poff+5]
        y[yoff+6] = p[poff+6]
    end

    # save aerodynamic inputs and distributed loads
    for i = 1:N
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_models = (section_aero, section_stru)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_stru)
        Nysi = number_of_inputs(section_stru)
        Npsi = number_of_parameters(section_stru)
        section_models = (section_aero, section_stru)

        # section properties
        icol = stru.icol_beam[i]
        element = assembly.elements[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab # local to body transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1] # local structural to aerodynamic transformation
        vi = R*CtCab'*Vinf - R*GXBeam.element_linear_velocity(element, P_elem, H_elem) # local linear freestream velocity
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem) # local angular freestream velocity

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        # calculate aerodynamic inputs and section loads

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ui = vcat(uai, usi)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pi = vcat(pai, psi)

        # section inputs
        yi = get_inputs(section_models, ui, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs

        # section aerodynamic loads (in body frame)
        fi = CtCab*R'*SVector(yi[Nyai+1], yi[Nyai+2], yi[Nyai+3])
        mi = CtCab*R'*SVector(yi[Nyai+4], yi[Nyai+5], yi[Nyai+6])

        # add constant distributed loads (in body frame)
        poff = 4 + 6*npoint + 6*(i-1)
        fi += SVector(padd[poff+1], padd[poff+2], padd[poff+3])
        mi += SVector(padd[poff+4], padd[poff+5], padd[poff+6])

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1 : yoff+3]] = fi
        y[iys[yoff+4 : yoff+6]] = mi
    end

    # save body frame linear/angular velocities
    yoff = nya + 6*npoint + 6*nelem
    y[yoff+1] = ur
    y[yoff+2] = vr
    y[yoff+3] = wr
    y[yoff+4] = pr
    y[yoff+5] = qr
    y[yoff+6] = rr

    return y
end

function get_input_mass_matrix!(My, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # zero out mass matrix
    My .= 0

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npadd = number_of_parameters(aero, stru) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    ipadd = npa + nps + 1 : npa + nps + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    us = view(u, ius) # structural state variables
    ps = view(p, ips) # structural parameters
    padd = view(p, ipadd) # additional parameters
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # global parameters
    Vx = padd[1]
    Vy = padd[2]
    Vz = padd[3]
    ρ = padd[4]
    ur = padd[end-5]
    vr = padd[end-4]
    wr = padd[end-3]
    pr = padd[end-2]
    qr = padd[end-1]
    rr = padd[end]

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - SVector(ur, vr, wr)

    # body angular velocity
    Ω = SVector(pr, qr, rr)

    # construct assembly from parameters
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # no point loads are dependent on state rates

    # save aerodynamic inputs and distributed loads
    for i = 1:N
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_models = (section_aero, section_stru)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_stru)
        Nysi = number_of_inputs(section_stru)
        Npsi = number_of_parameters(section_stru)
        section_models = (section_aero, section_stru)

        # section properties
        icol = stru.icol_beam[i]
        element = assembly.elements[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab # local to body transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1] # local structural to aerodynamic transformation
        vi = R*CtCab'*Vinf - R*GXBeam.element_linear_velocity(element, P_elem, H_elem) # local linear freestream velocity
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem) # local angular freestream velocity
        dvi_dPi = -R * element.minv11 * stru.mass_scaling
        dvi_dHi = -R * element.minv12 * stru.mass_scaling
        dωi_dPi = R * element.minv12' * stru.mass_scaling
        dωi_dHi = R * element.minv22 * stru.mass_scaling

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up

        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ui = vcat(uai, usi)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pi = vcat(pai, psi)

        # section input mass matrix
        Myi = get_input_mass_matrix(section_models, ui, pi, t)

        # separate into component mass matrices
        yai_duai = SMatrix{Nyai,Nuai}(view(Myi, 1:Nyai, 1:Nuai))
        yai_dvi = SMatrix{Nyai,3}(view(Myi, 1:Nyai, Nuai+1:Nuai+3))
        yai_dωi = SMatrix{Nyai,3}(view(Myi, 1:Nyai, Nuai+4:Nuai+6))
        f_duai = SMatrix{3,Nuai}(view(Myi, Nyai+1:Nyai+3, 1:Nuai))
        f_dvi = SMatrix{3,3}(view(Myi, Nyai+1:Nyai+3, Nuai+1:Nuai+3))
        f_dωi = SMatrix{3,3}(view(Myi, Nyai+1:Nyai+3, Nuai+4:Nuai+6))
        m_duai = SMatrix{3,Nuai}(view(Myi, Nyai+4:Nyai+6, 1:Nuai))
        m_dvi = SMatrix{3,3}(view(Myi, Nyai+4:Nyai+6, Nuai+1:Nuai+3))
        m_dωi = SMatrix{3,3}(view(Myi, Nyai+4:Nyai+6, Nuai+4:Nuai+6))

        # propagate derivatives using chain rule
        yai_dPi = yai_dvi * dvi_dPi + yai_dωi * dωi_dPi
        yai_dHi = yai_dvi * dvi_dHi + yai_dωi * dωi_dHi
        f_dPi = f_dvi * dvi_dPi + f_dωi * dωi_dPi
        f_dHi = f_dvi * dvi_dHi + f_dωi * dωi_dHi
        m_dPi = m_dvi * dvi_dPi + m_dωi * dωi_dPi
        m_dHi = m_dvi * dvi_dHi + m_dωi * dωi_dHi

        # save local inputs
        icol = stru.icol_beam[i]

        # aerodynamic inputs
        My[iya[iyas[i]], iua[iuas[i]]] = yai_duai
        My[iya[iyas[i]], ius[icol+12]:ius[icol+14]] = yai_dPi
        My[iya[iyas[i]], ius[icol+15]:ius[icol+17]] = yai_dHi

        # save load mass matrix entries (in body frame)
        offset = 6*npoint + 6*(i-1)
        My[iys[offset+1 : offset+3], iua[iuas[i]]] = CtCab*R'*f_duai
        My[iys[offset+4 : offset+6], iua[iuas[i]]] = CtCab*R'*m_duai
        My[iys[offset+1 : offset+3], ius[icol+12]:ius[icol+14]] = CtCab*R'*f_dPi
        My[iys[offset+1 : offset+3], ius[icol+15]:ius[icol+17]] = CtCab*R'*f_dHi
        My[iys[offset+4 : offset+6], ius[icol+12]:ius[icol+14]] = CtCab*R'*m_dPi
        My[iys[offset+4 : offset+6], ius[icol+15]:ius[icol+17]] = CtCab*R'*m_dHi
    end

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{N,T}, stru::GEBT,
    du, u, p, t) where {N,T}

    # initialize input vector
    models = (aero, stru)
    TF = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    y = zeros(TF, Ny)

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npadd = number_of_parameters(aero, stru) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    ipadd = npa + nps + 1 : npa + nps + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state rates, states, inputs, and parameters
    dua = view(du, iua) # aerodynamic states rates
    ua = view(u, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dus = view(du, ius) # structural state rates
    us = view(u, ius) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic state rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # global parameters
    Vx = padd[1]
    Vy = padd[2]
    Vz = padd[3]
    ρ = padd[4]
    ur = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 1]
    vr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 2]
    wr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 3]
    pr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 4]
    qr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 5]
    rr = padd[4 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam) + 6]

    # rigid body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - V

    # construct assembly from parameters
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # no point loads are dependent on state rates

    # save aerodynamic inputs and distributed loads
    for i = 1:N
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_models = (section_aero, section_stru)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_aero)
        Nysi = number_of_inputs(section_aero)
        Npsi = number_of_parameters(section_aero)

        # local section properties
        icol = stru.icol_beam[i]
        element = assembly.elements[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum

        du_elem = SVector(dus[icol], dus[icol+1], dus[icol+2]) # linear displacement
        dθ_elem = SVector(dus[icol+3], dus[icol+4], dus[icol+5]) # angular displacement
        dF_elem = SVector(dus[icol+6], dus[icol+7], dus[icol+8]) .* stru.force_scaling # internal force
        dM_elem = SVector(dus[icol+9], dus[icol+10], dus[icol+11]) .* stru.force_scaling # internal moment
        dP_elem = SVector(dus[icol+12], dus[icol+13], dus[icol+14]) .* stru.mass_scaling # linear momentum
        dH_elem = SVector(dus[icol+15], dus[icol+16], dus[icol+17]) .* stru.mass_scaling # angular momentum

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)
        dθ_elem *= scaling # angular displacement rate (Wiener-Milenkovic parameters)

        # local to body transformation
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab

        # local structural to aerodynamic transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # local freestream linear and angular velocities
        vi = R*CtCab'*Vinf - R*GXBeam.element_linear_velocity(element, P_elem, H_elem)
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        # local freestream linear and angular accelerations
        dvi = -R*GXBeam.element_linear_velocity(element, dP_elem, dH_elem)
        dωi = R*GXBeam.element_angular_velocity(element, dP_elem, dH_elem)

        # section state variable rates
        duai = SVector{Nuai}(duas[i]) # aerodynamic state variables
        dusi = vcat(dvi, dωi) # structural state variables
        dui = vcat(duai, dusi)

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ui = vcat(uai, usi)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pi = vcat(pai, psi)

        # section inputs (from state rates)
        yi = get_inputs_from_state_rates(section_models, dui, ui, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs

        # section aerodynamic loads (in body frame)
        fi = CtCab*R'*SVector(yi[Nyai+1], yi[Nyai+2], yi[Nyai+3])
        mi = CtCab*R'*SVector(yi[Nyai+4], yi[Nyai+5], yi[Nyai+6])

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1 : yoff+3]] = fi
        y[iys[yoff+4 : yoff+6]] = mi

    end

    return y
end
