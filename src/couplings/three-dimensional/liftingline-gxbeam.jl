"""
    couple_models(aero::LiftingLine, stru::GEBT)

Create a coupled model using a lifting line aerodynamic model and a
geometrically exact beam theory model.  This model introduces additional
parameters corresponding to the freestream air density ``\\rho_\\infty``,
followed by the external loads ``F_{x,i}, F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i},
M_{z,i}`` or displacements ``u_{x,i}, u_{y,i}, u_{z,i}, \\theta_{x,i},
\\theta_{y,i}, \\theta_{z,i}`` for each node, followed by the constant
distributed loads ``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i}, m_{y,i}, m_{z,i}``
applied on each beam element (excluding aerodynamic loads), followed by the
body frame linear and angular velocities ``u, v, w, p, q, r``.

**NOTE: When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface**
"""
couple_models(aero::LiftingLine, stru::GEBT) = (aero, stru)

# --- traits --- #

function number_of_additional_parameters(aero::LiftingLine, stru::GEBT)
    return 1 + 6*length(stru.icol_point) + 6*length(stru.icol_elem) + 6
end

function coupling_inplaceness(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT}) where {NA,TA}
    return InPlace()
end

function coupling_rate_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT}) where {NA,TA}
    model_types = TA.parameters
    if all(isempty.(coupling_rate_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(coupling_rate_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(coupling_rate_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isinvariant.(coupling_rate_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Invariant()
    elseif all(isconstant.(coupling_rate_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(coupling_rate_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function coupling_state_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT}) where {NA,TA}
    model_types = TA.parameters
    if all(isempty.(coupling_state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(coupling_state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isidentity.(coupling_state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Identity()
    elseif all(isinvariant.(coupling_state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Invariant()
    elseif all(isconstant.(coupling_state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(coupling_state_jacobian_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

coupling_parameter_jacobian_type(::Type{<:LiftingLine}, ::Type{<:GEBT}) = Nonlinear()

function coupling_time_gradient_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT}) where {NA,TA}
    model_types = TA.parameters
    if all(isempty.(coupling_time_gradient_type.(model_types, Ref(LiftingLineSection))))
        return Empty()
    elseif all(iszero.(coupling_time_gradient_type.(model_types, Ref(LiftingLineSection))))
        return Zeros()
    elseif all(isinvariant.(coupling_time_gradient_type.(model_types, Ref(LiftingLineSection))))
        return Invariant()
    elseif all(isconstant.(coupling_time_gradient_type.(model_types, Ref(LiftingLineSection))))
        return Constant()
    elseif all(islinear.(coupling_time_gradient_type.(model_types, Ref(LiftingLineSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

# --- methods --- #

function get_coupling_inputs!(y, aero::LiftingLine{NA,TA}, stru::GEBT, dx, x, p, t) where {NA,TA}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npadd = number_of_additional_parameters(aero, stru) # number of additional parameters

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
    dua = view(dx, iua) # aerodynamic rates
    ua = view(x, iua) # aerodynamic states
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dus = view(dx, ius) # structural rates
    us = view(x, ius) # structural states
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic states for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # number of points and elements
    npoint = length(stru.icol_point)
    nelem = length(stru.icol_elem)

    # global parameters
    ρ = padd[1]
    ur = padd[1 + 6*npoint + 6*nelem + 1]
    vr = padd[1 + 6*npoint + 6*nelem + 2]
    wr = padd[1 + 6*npoint + 6*nelem + 3]
    pr = padd[1 + 6*npoint + 6*nelem + 4]
    qr = padd[1 + 6*npoint + 6*nelem + 5]
    rr = padd[1 + 6*npoint + 6*nelem + 6]

    # body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # construct assembly from structural parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # set prescribed point loads/displacements
    for ip = 1:npoint
        yoff = 6*(ip-1)
        poff = 1 + 6*(ip-1)
        ys[yoff+1] = padd[poff+1]
        ys[yoff+2] = padd[poff+2]
        ys[yoff+3] = padd[poff+3]
        ys[yoff+4] = padd[poff+4]
        ys[yoff+5] = padd[poff+5]
        ys[yoff+6] = padd[poff+6]
    end

    # set aerodynamic inputs and distributed loads
    for i = 1:NA
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

        # section properties
        icol = stru.icol_elem[i]
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
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab # local to body transformation

        # local structural to aerodynamic transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1] # local structural to aerodynamic transformation

        # local freestream linear and angular velocities
        vi = -R*(CtCab'*V + GXBeam.element_linear_velocity(element, P_elem, H_elem)) # local linear freestream velocity
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem) # local angular freestream velocity

        # local freestream linear and angular accelerations
        dvi = -R*GXBeam.element_linear_velocity(element, dP_elem, dH_elem)
        dωi = R*GXBeam.element_angular_velocity(element, dP_elem, dH_elem)

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

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

        # section inputs
        yi = get_coupling_inputs(section_models, dui, ui, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs

        # section aerodynamic loads (in body frame)
        fi = CtCab*R'*SVector(ysi[1], ysi[2], ysi[3])
        mi = CtCab*R'*SVector(ysi[4], ysi[5], ysi[6])

        # add constant distributed loads (in body frame)
        poff = 1 + 6*npoint + 6*(i-1)
        fi += SVector(padd[poff+1], padd[poff+2], padd[poff+3])
        mi += SVector(padd[poff+4], padd[poff+5], padd[poff+6])

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1 : yoff+3]] = fi
        y[iys[yoff+4 : yoff+6]] = mi
    end

    # save body frame linear/angular velocities
    yoff = 6*npoint + 6*nelem
    y[iys[yoff+1]] = ur
    y[iys[yoff+2]] = vr
    y[iys[yoff+3]] = wr
    y[iys[yoff+4]] = pr
    y[iys[yoff+5]] = qr
    y[iys[yoff+6]] = rr

    return y
end

# --- performance overloads --- #

# TODO

# --- convenience methods --- #

function set_additional_parameters!(padd, aero::LiftingLine, stru::GEBT;
    rho, point_conditions, element_loads, u, v, w, p, q, r) where {NA,TA}

    np = length(stru.icol_point)
    ne = length(stru.icol_elem)

    padd[1] = rho

    for ip = 1:np
        padd[1+6*(ip-1)+1] = point_conditions[1,ip]
        padd[1+6*(ip-1)+2] = point_conditions[2,ip]
        padd[1+6*(ip-1)+3] = point_conditions[3,ip]
        padd[1+6*(ip-1)+4] = point_conditions[4,ip]
        padd[1+6*(ip-1)+5] = point_conditions[5,ip]
        padd[1+6*(ip-1)+6] = point_conditions[6,ip]
    end

    for ie = 1:ne
        padd[1+6*np+6*(ie-1)+1] = element_loads[1,ie]
        padd[1+6*np+6*(ie-1)+2] = element_loads[2,ie]
        padd[1+6*np+6*(ie-1)+3] = element_loads[3,ie]
        padd[1+6*np+6*(ie-1)+4] = element_loads[4,ie]
        padd[1+6*np+6*(ie-1)+5] = element_loads[5,ie]
        padd[1+6*np+6*(ie-1)+6] = element_loads[6,ie]
    end

    padd[1+6*np+6*ne+1] = u
    padd[1+6*np+6*ne+2] = v
    padd[1+6*np+6*ne+3] = w
    padd[1+6*np+6*ne+4] = p
    padd[1+6*np+6*ne+5] = q
    padd[1+6*np+6*ne+6] = r

    return padd
end

function separate_additional_parameters(aero::LiftingLine, stru::GEBT, padd)

    np = length(stru.icol_point)
    ne = length(stru.icol_elem)

    rho = padd[1]

    point_conditions = reshape(view(padd, 1 + 1 : 1 + 6*np), 6, np)

    element_loads = reshape(view(padd, 1 + 6*np + 1 : 1 + 6*np + 6*ne), 6, ne)

    u = padd[1 + 6*np + 6*ne + 1]
    v = padd[1 + 6*np + 6*ne + 2]
    w = padd[1 + 6*np + 6*ne + 3]
    p = padd[1 + 6*np + 6*ne + 4]
    q = padd[1 + 6*np + 6*ne + 5]
    r = padd[1 + 6*np + 6*ne + 6]

    return (rho = rho, point_conditions = point_conditions, element_loads = element_loads,
        u = u, v = v, w = w, p = p, q = q, r = r)
end
