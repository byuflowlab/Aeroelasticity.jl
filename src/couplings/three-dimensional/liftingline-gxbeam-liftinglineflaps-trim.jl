"""
    couple_models(aero::LiftingLine, stru::GEBT, flap::LiftingLineFlaps, ctrl::Trim)

Create a coupled model using a lifting line aerodynamic model, a
geometrically exact beam theory model, a lifting line control surface model,
and a controller which maintains trimmed operating conditions.  This model introduces
additional parameters corresponding to the freestream air density ``\\rho_\\infty``,
followed by the external loads ``F_{x,i}, F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i},
M_{z,i}`` or displacements ``u_{x,i}, u_{y,i}, u_{z,i}, \\theta_{x,i},
\\theta_{y,i}, \\theta_{z,i}`` for each node, followed by the constant distributed
loads ``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i}, m_{y,i}, m_{z,i}`` for each beam
element (excluding aerodynamic loads), followed by the body frame linear and
angular velocities ``u, v, w, p, q, r``, and the control surface deflections
which are not defined by the controller.

**NOTE: When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface**
"""
function couple_models(aero::LiftingLine, stru::GEBT, flap::LiftingLineFlaps,
    ctrl::Trim)

    return (aero, stru, flap, ctrl)
end

# --- traits --- #


function number_of_additional_parameters(aero::LiftingLine, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}) where {NA,NF,NG,NC,TF,TG}

    return 1 + 6*length(stru.icol_point) + 6*length(stru.icol_elem) + 6 + NG - NC
end

function coupling_inplaceness(::Type{<:LiftingLine}, ::Type{<:GEBT},
    ::Type{<:LiftingLineFlaps}, ::Type{<:Trim})

    return InPlace()
end

function coupling_mass_matrix_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}, ::Type{<:Trim}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(isempty.(coupling_mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Empty()
    elseif all(iszero.(coupling_mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Zeros()
    elseif all(isidentity.(coupling_mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Identity()
    elseif all(isconstant.(coupling_mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Constant()
    elseif all(islinear.(coupling_mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function coupling_state_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}, ::Type{<:Trim}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(isempty.(coupling_state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Empty()
    elseif all(iszero.(coupling_state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Zeros()
    elseif all(isidentity.(coupling_state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Identity()
    elseif all(isconstant.(coupling_state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Constant()
    elseif all(islinear.(coupling_state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{NA,TA}, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, u, p, t) where {NA,NF,NG,NC,TA,TF,TG}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nuf = number_of_states(flap) # number of control surface states
    nuc = number_of_states(ctrl) # number of control surface states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyf = number_of_inputs(flap) # number of control surface inputs
    nyc = number_of_inputs(ctrl) # number of control surface inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npf = number_of_parameters(flap) # number of control surface parameters
    npc = number_of_parameters(ctrl) # number of control surface parameters
    npadd = number_of_additional_parameters(aero, stru, flap, ctrl) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iuf = nua + nus + 1 : nua + nus + nuf # indices of control surface states
    iyf = nya + nys + 1 : nya + nys + nyf # indices of control surface inputs
    ipf = npa + nps + 1 : npa + nps + npf # indices of control surface parameters
    iuc = nua + nus + nuf + 1 : nua + nus + nuf + nuc # indices of control surface states
    iyc = nya + nys + nyf + 1 : nya + nys + nyf + nyc # indices of control surface inputs
    ipc = npa + nps + npf + 1 : npa + nps + npf + npc # indices of control surface parameters
    ipadd = npa + nps + npf + npc + 1 : npa + nps + npf + npc + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic parameters for each section
    iufs = state_indices(flap.models) # indices of control surface states for each section
    iyfs = input_indices(flap.models) # indices of control surface inputs for each section
    ipfs = parameter_indices(flap.models) # indices of control surface parameters for each section

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    us = view(u, ius) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    uf = view(u, iuf) # structural state variables
    yf = view(y, iyf) # structural inputs
    pf = view(p, ipf) # structural parameters
    uc = view(u, iuc) # structural state variables
    yc = view(y, iyc) # structural inputs
    pc = view(p, ipc) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # number of points and elements
    npoint = length(stru.icol_point)
    nelem = length(stru.icol_elem)

    # controller state variables
    δc = uc

    # global parameters
    ρ = padd[1]
    ur = padd[1 + 6*npoint + 6*nelem + 1]
    vr = padd[1 + 6*npoint + 6*nelem + 2]
    wr = padd[1 + 6*npoint + 6*nelem + 3]
    pr = padd[1 + 6*npoint + 6*nelem + 4]
    qr = padd[1 + 6*npoint + 6*nelem + 5]
    rr = padd[1 + 6*npoint + 6*nelem + 6]
    δp = view(padd, 1 + 6*npoint + 6*nelem + 6 + 1 : 1 + 6*npoint + 6*nelem + 6 + NG - NC)

    # commanded control deflections
    T = promote_type(eltype(u), eltype(p))
    δ = zeros(SVector{NG,T})
    iδu = iδp = 1
    for i = 1:NG
        if i in ctrl.state_indices
            δ = setindex(δ, δc[iδu], i)
            iδu += 1
        else
            δ = setindex(δ, δp[iδp], i)
            iδp += 1
        end
    end

    # body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # construct assembly from structural parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # save prescribed point loads/displacements
    for ip = 1:npoint
        yoff = 6*(ip-1)
        poff = 1 + 6*(ip-1)
        y[iys[yoff+1]] = padd[poff+1]
        y[iys[yoff+2]] = padd[poff+2]
        y[iys[yoff+3]] = padd[poff+3]
        y[iys[yoff+4]] = padd[poff+4]
        y[iys[yoff+5]] = padd[poff+5]
        y[iys[yoff+6]] = padd[poff+6]
    end

    # add prescribed loads to total forces and moments
    for ip = 1:npoint
        # point location
        point = assembly.points[ip]
        # prescribed condition identities
        prescribed_forces = SVector{6}(view(stru.displacement, :, ip)) .== false
        # point loads
        poff = 1 + 6*(ip-1)
        F1 = ifelse(prescribed_forces[1], padd[poff+1], zero(padd[poff+1]))
        F2 = ifelse(prescribed_forces[2], padd[poff+2], zero(padd[poff+2]))
        F3 = ifelse(prescribed_forces[3], padd[poff+3], zero(padd[poff+3]))
        M1 = ifelse(prescribed_forces[4], padd[poff+4], zero(padd[poff+4]))
        M2 = ifelse(prescribed_forces[5], padd[poff+5], zero(padd[poff+5]))
        M3 = ifelse(prescribed_forces[6], padd[poff+6], zero(padd[poff+6]))
        Fp = SVector(F1, F2, F3)
        Mp = SVector(M1, M2, M3)
        Ftot += Fp
        Mtot += cross(point, Fp) + Mp
    end

    # calculate aerodynamic inputs, distributed loads, and element properties
    for i = 1:NA
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_flap = flap.models[i]
        section_ctrl = LiftingLineSectionControl()
        section_models = (section_aero, section_stru, section_flap, section_ctrl)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_stru)
        Nysi = number_of_inputs(section_stru)
        Npsi = number_of_parameters(section_stru)
        Nufi = number_of_states(section_flap)
        Nyfi = number_of_inputs(section_flap)
        Npfi = number_of_parameters(section_flap)
        Nuci = number_of_states(section_ctrl)
        Nyci = number_of_inputs(section_ctrl)
        Npci = number_of_parameters(section_ctrl)

        # local section properties
        icol = stru.icol_elem[i]
        element = assembly.elements[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)
        ΔL = element.L # element length
        pe = element.x + u_elem # location of element (in body frame)
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab # local to body transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1] # local structural to aerodynamic transformation
        vi = -R*(CtCab'*V + R*GXBeam.element_linear_velocity(element, P_elem, H_elem)) # local linear freestream velocity
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem) # local angular freestream velocity

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = SMatrix{1,NG}(getindex.(flap.gains, i)) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pfi = SVector{Npfi}(pfs[i]) # control surface parameters
        pci = SVector{0,Float64}() # controller parameters
        pi = vcat(pai, psi, pfi, pci)

        # section inputs
        yi = get_inputs(section_models, ui, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)
        yfi = view(yi, Nyai + Nysi + 1 : Nyai + Nysi + Nyfi)
        yci = view(yi, Nyai + Nysi + Nyfi + 1 : Nyai + Nysi + Nyfi + Nyci)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs
        y[iyf[iyfs[i]]] .= yfi # control surface inputs

        # section aerodynamic loads (in body frame)
        fi = CtCab*R'*SVector(yi[Nyai+1], yi[Nyai+2], yi[Nyai+3])
        mi = CtCab*R'*SVector(yi[Nyai+4], yi[Nyai+5], yi[Nyai+6])

        # add constant distributed loads (in body frame)
        poff = 1 + 6*npoint + 6*(i-1)
        fi += SVector(padd[poff+1], padd[poff+2], padd[poff+3])
        mi += SVector(padd[poff+4], padd[poff+5], padd[poff+6])

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1 : yoff+3]] = fi
        y[iys[yoff+4 : yoff+6]] = mi

        # add distributed loads to total forces and moments
        Ftot += ΔL*fi
        Mtot += cross(pe, ΔL*fi) + ΔL*mi

    end

    # save body frame linear/angular velocities
    yoff = 6*npoint + 6*nelem
    y[iys[yoff+1]] = ur
    y[iys[yoff+2]] = vr
    y[iys[yoff+3]] = wr
    y[iys[yoff+4]] = pr
    y[iys[yoff+5]] = qr
    y[iys[yoff+6]] = rr

    # save trim model inputs
    loads = SVector(Ftot[1], Ftot[2], Ftot[3], Mtot[1], Mtot[2], Mtot[3])

    for i = 1:NC
        idx = ctrl.force_indices[i]
        y[iyc[i]] = loads[idx]
    end

    return y
end

function get_coupling_mass_matrix!(My, aero::LiftingLine{NA,TA}, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, u, p, t) where {NA,NF,NG,NC,TA,TF,TG}

    # zero out mass matrix
    My .= 0

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nuf = number_of_states(flap) # number of control surface states
    nuc = number_of_states(ctrl) # number of control surface states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyf = number_of_inputs(flap) # number of control surface inputs
    nyc = number_of_inputs(ctrl) # number of control surface inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npf = number_of_parameters(flap) # number of control surface parameters
    npc = number_of_parameters(ctrl) # number of control surface parameters
    npadd = number_of_additional_parameters(aero, stru, flap, ctrl) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iuf = nua + nus + 1 : nua + nus + nuf # indices of control surface states
    iyf = nya + nys + 1 : nya + nys + nyf # indices of control surface inputs
    ipf = npa + nps + 1 : npa + nps + npf # indices of control surface parameters
    iuc = nua + nus + nuf + 1 : nua + nus + nuf + nuc # indices of control surface states
    iyc = nya + nys + nyf + 1 : nya + nys + nyf + nyc # indices of control surface inputs
    ipc = npa + nps + npf + 1 : npa + nps + npf + npc # indices of control surface parameters
    ipadd = npa + nps + npf + npc + 1 : npa + nps + npf + npc + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic parameters for each section
    iufs = state_indices(flap.models) # indices of control surface states for each section
    iyfs = input_indices(flap.models) # indices of control surface inputs for each section
    ipfs = parameter_indices(flap.models) # indices of control surface parameters for each section

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    us = view(u, ius) # structural state variables
    ps = view(p, ips) # structural parameters
    uf = view(u, iuf) # structural state variables
    pf = view(p, ipf) # structural parameters
    uc = view(u, iuc) # structural state variables
    pc = view(p, ipc) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # number of points and elements
    npoint = length(stru.icol_point)
    nelem = length(stru.icol_elem)

    # control state variables
    δc = uc

    # global parameters
    ρ = padd[1]
    ur = padd[1 + 6*npoint + 6*nelem + 1]
    vr = padd[1 + 6*npoint + 6*nelem + 2]
    wr = padd[1 + 6*npoint + 6*nelem + 3]
    pr = padd[1 + 6*npoint + 6*nelem + 4]
    qr = padd[1 + 6*npoint + 6*nelem + 5]
    rr = padd[1 + 6*npoint + 6*nelem + 6]
    δp = view(padd, 1 + 6*npoint + 6*nelem + 6 + 1 : 1 + 6*npoint + 6*nelem + 6 + NG - NC)

    # commanded control deflections
    T = promote_type(eltype(u), eltype(p))
    δ = zeros(SVector{NG,T})
    iδu = iδp = 1
    for i = 1:NG
        if i in ctrl.state_indices
            δ = setindex(δ, δc[iδu], i)
            iδu += 1
        else
            δ = setindex(δ, δp[iδp], i)
            iδp += 1
        end
    end

    # rigid body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # no point loads are dependent on state rates

    # calculate aerodynamic inputs, distributed loads, and element properties
    for i = 1:NA
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_flap = flap.models[i]
        section_ctrl = LiftingLineSectionControl()
        section_models = (section_aero, section_stru, section_flap, section_ctrl)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_stru)
        Nysi = number_of_inputs(section_stru)
        Npsi = number_of_parameters(section_stru)
        Nufi = number_of_states(section_flap)
        Nyfi = number_of_inputs(section_flap)
        Npfi = number_of_parameters(section_flap)
        Nuci = number_of_states(section_ctrl)
        Nyci = number_of_inputs(section_ctrl)
        Npci = number_of_parameters(section_ctrl)

        # local section properties
        icol = stru.icol_elem[i]
        element = assembly.elements[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)
        ΔL = element.L # element length
        pe = element.x + u_elem # location of element (in body frame)
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab # local to body transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1] # local structural to aerodynamic transformation
        vi = -R*(CtCab'*V + GXBeam.element_linear_velocity(element, P_elem, H_elem)) # local linear freestream velocity
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
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = SMatrix{1,NG}(getindex.(flap.gains, i)) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pfi = SVector{Npfi}(pfs[i]) # control surface parameters
        pci = SVector{0,Float64}() # controller parameters
        pi = vcat(pai, psi, pfi, pci)

        # section input mass matrix
        Myi = get_coupling_mass_matrix(section_models, ui, pi, t)

        # separate into component mass matrices
        yai_duai = SMatrix{Nyai,Nuai}(view(Myi, 1:Nyai, 1:Nuai))
        yai_dvi = SMatrix{Nyai,3}(view(Myi, 1:Nyai, Nuai+1:Nuai+3))
        yai_dωi = SMatrix{Nyai,3}(view(Myi, 1:Nyai, Nuai+4:Nuai+6))
        yai_dufi = SMatrix{Nyai,Nufi}(view(Myi, 1:Nyai, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        yai_duci = SMatrix{Nyai,Nuci}(view(Myi, 1:Nyai, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        f_duai = SMatrix{3,Nuai}(view(Myi, Nyai+1:Nyai+3, 1:Nuai))
        f_dvi = SMatrix{3,3}(view(Myi, Nyai+1:Nyai+3, Nuai+1:Nuai+3))
        f_dωi = SMatrix{3,3}(view(Myi, Nyai+1:Nyai+3, Nuai+4:Nuai+6))
        f_dufi = SMatrix{3,Nufi}(view(Myi, Nyai+1:Nyai+3, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        f_duci = SMatrix{3,Nuci}(view(Myi, Nyai+1:Nyai+3, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        m_duai = SMatrix{3,Nuai}(view(Myi, Nyai+4:Nyai+6, 1:Nuai))
        m_dvi = SMatrix{3,3}(view(Myi, Nyai+4:Nyai+6, Nuai+1:Nuai+3))
        m_dωi = SMatrix{3,3}(view(Myi, Nyai+4:Nyai+6, Nuai+4:Nuai+6))
        m_dufi = SMatrix{3,Nufi}(view(Myi, Nyai+4:Nyai+6, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        m_duci = SMatrix{3,Nuci}(view(Myi, Nyai+4:Nyai+6, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        yfi_duai = SMatrix{Nyfi,Nuai}(view(Myi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, 1:Nuai))
        yfi_dvi = SMatrix{Nyfi,3}(view(Myi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+1:Nuai+3))
        yfi_dωi = SMatrix{Nyfi,3}(view(Myi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+4:Nuai+6))
        yfi_dufi = SMatrix{Nyfi,Nufi}(view(Myi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        yfi_duci = SMatrix{Nyfi,Nuci}(view(Myi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        yci_duai = SMatrix{Nyci,Nuai}(view(Myi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, 1:Nuai))
        yci_dvi = SMatrix{Nyci,3}(view(Myi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+1:Nuai+3))
        yci_dωi = SMatrix{Nyci,3}(view(Myi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+4:Nuai+6))
        yci_dufi = SMatrix{Nyci,Nufi}(view(Myi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        yci_duci = SMatrix{Nyci,Nuci}(view(Myi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))

        # propagate derivatives using chain rule
        yai_dPi = yai_dvi * dvi_dPi + yai_dωi * dωi_dPi
        yai_dHi = yai_dvi * dvi_dHi + yai_dωi * dωi_dHi
        f_dPi = f_dvi * dvi_dPi + f_dωi * dωi_dPi
        f_dHi = f_dvi * dvi_dHi + f_dωi * dωi_dHi
        m_dPi = m_dvi * dvi_dPi + m_dωi * dωi_dPi
        m_dHi = m_dvi * dvi_dHi + m_dωi * dωi_dHi
        yfi_dPi = yfi_dvi * dvi_dPi + yfi_dωi * dωi_dPi
        yfi_dHi = yfi_dvi * dvi_dHi + yfi_dωi * dωi_dHi

        # save local inputs
        icol = stru.icol_elem[i]

        # aerodynamic inputs
        My[iya[iyas[i]], iua[iuas[i]]] = yai_duai
        My[iya[iyas[i]], ius[icol+12]:ius[icol+14]] = yai_dPi
        My[iya[iyas[i]], ius[icol+15]:ius[icol+17]] = yai_dHi
        My[iya[iyas[i]], iuf[iufs[i]]] = yai_dufi

        # control surface inputs
        My[iyf[iyfs[i]], iua[iuas[i]]] = yfi_duai
        My[iyf[iyfs[i]], ius[icol+12]:ius[icol+14]] = yfi_dPi
        My[iyf[iyfs[i]], ius[icol+15]:ius[icol+17]] = yfi_dHi
        My[iyf[iyfs[i]], iuf[iufs[i]]] = yfi_dufi

        # section aerodynamic loads (in body frame)
        fi_duai = CtCab*R'*f_duai
        fi_dPi = CtCab*R'*f_dPi
        fi_dHi = CtCab*R'*f_dHi
        fi_dufi = CtCab*R'*f_dufi
        mi_duai = CtCab*R'*m_duai
        mi_dPi = CtCab*R'*m_dPi
        mi_dHi = CtCab*R'*m_dHi
        mi_dufi = CtCab*R'*m_dufi

        # save load mass matrix entries (in body frame)
        yoff = 6*npoint + 6*(i-1)
        My[iys[yoff+1 : yoff+3], iua[iuas[i]]] = fi_duai
        My[iys[yoff+4 : yoff+6], iua[iuas[i]]] = mi_duai
        My[iys[yoff+1 : yoff+3], ius[icol+12]:ius[icol+14]] = fi_dPi
        My[iys[yoff+1 : yoff+3], ius[icol+15]:ius[icol+17]] = fi_dHi
        My[iys[yoff+4 : yoff+6], ius[icol+12]:ius[icol+14]] = mi_dPi
        My[iys[yoff+4 : yoff+6], ius[icol+15]:ius[icol+17]] = mi_dHi
        My[iys[yoff+1 : yoff+3], iuf[iufs[i]]] = fi_dufi
        My[iys[yoff+4 : yoff+6], iuf[iufs[i]]] = mi_dufi

        loads_duai = vcat(ΔL*fi_duai, GXBeam.tilde(pe)*ΔL*fi_duai + ΔL*mi_duai)
        loads_dPi = vcat(ΔL*fi_dPi, GXBeam.tilde(pe)*ΔL*fi_dPi + ΔL*mi_dPi)
        loads_dHi = vcat(ΔL*fi_dHi, GXBeam.tilde(pe)*ΔL*fi_dHi + ΔL*mi_dHi)
        loads_dufi = vcat(ΔL*fi_dufi, GXBeam.tilde(pe)*ΔL*fi_dufi + ΔL*mi_dufi)

        for ic = 1:NC
            idx = ctrl.force_indices[ic]
            My[iyc[ic], iua[iuas[i]]] = loads_duai[idx,:]
            My[iyc[ic], ius[icol+12]:ius[icol+14]] = loads_dPi[idx,:]
            My[iyc[ic], ius[icol+15]:ius[icol+17]] = loads_dHi[idx,:]
            My[iyc[ic], iuf[iufs[i]]] = loads_dufi[idx,:]
        end
    end

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_using_state_rates(aero::LiftingLine{NA,TA}, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, du, u, p, t) where {NA,NF,NG,NC,TA,TF,TG}

    # initialize input vector
    models = (aero, stru, flap, ctrl)
    T = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    y = zeros(T, Ny)

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nuf = number_of_states(flap) # number of control surface states
    nuc = number_of_states(ctrl) # number of control surface states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyf = number_of_inputs(flap) # number of control surface inputs
    nyc = number_of_inputs(ctrl) # number of control surface inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npf = number_of_parameters(flap) # number of control surface parameters
    npc = number_of_parameters(ctrl) # number of control surface parameters
    npadd = number_of_additional_parameters(aero, stru, flap, ctrl) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iuf = nua + nus + 1 : nua + nus + nuf # indices of control surface states
    iyf = nya + nys + 1 : nya + nys + nyf # indices of control surface inputs
    ipf = npa + nps + 1 : npa + nps + npf # indices of control surface parameters
    iuc = nua + nus + nuf + 1 : nua + nus + nuf + nuc # indices of control surface states
    iyc = nya + nys + nyf + 1 : nya + nys + nyf + nyc # indices of control surface inputs
    ipc = npa + nps + npf + 1 : npa + nps + npf + npc # indices of control surface parameters
    ipadd = npa + nps + npf + npc + 1 : npa + nps + npf + npc + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic parameters for each section
    iufs = state_indices(flap.models) # indices of control surface states for each section
    iyfs = input_indices(flap.models) # indices of control surface inputs for each section
    ipfs = parameter_indices(flap.models) # indices of control surface parameters for each section

    # separate state variables, inputs, and parameters
    dua = view(du, iua) # aerodynamic state variables
    ua = view(u, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dus = view(du, ius) # structural state variables
    us = view(u, ius) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    duf = view(du, iuf) # control surface state variables
    uf = view(u, iuf) # control surface state variables
    yf = view(y, iyf) # control surface inputs
    pf = view(p, ipf) # control surface parameters
    duc = view(du, iuc) # controller state rates
    uc = view(u, iuc) # controller state variables
    yc = view(y, iyc) # controller inputs
    pc = view(p, ipc) # controller parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic state variables for each section
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    dufs = view.(Ref(duf), iufs) # control surface state variables for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # number of points and elements
    npoint = length(stru.icol_point)
    nelem = length(stru.icol_elem)

    # control state variables
    δc = uc

    # global parameters
    ρ = padd[1]
    ur = padd[1 + 6*npoint + 6*nelem + 1]
    vr = padd[1 + 6*npoint + 6*nelem + 2]
    wr = padd[1 + 6*npoint + 6*nelem + 3]
    pr = padd[1 + 6*npoint + 6*nelem + 4]
    qr = padd[1 + 6*npoint + 6*nelem + 5]
    rr = padd[1 + 6*npoint + 6*nelem + 6]
    δp = view(padd, 1 + 6*npoint + 6*nelem + 6 + 1 : 1 + 6*npoint + 6*nelem + 6 + NG - NC)

    # commanded control deflections
    T = promote_type(eltype(du), eltype(u), eltype(p))
    dδ = zeros(SVector{NG,T})
    iδu = 1
    for i = 1:NG
        if i in ctrl.state_indices
            dδ = setindex(dδ, δc[iδu], i)
            iδu += 1
        else
            dδ = setindex(dδ, 0, i)
        end
    end

    # commanded control deflections
    T = promote_type(eltype(u), eltype(p))
    δ = zeros(SVector{NG,T})
    iδu = iδp = 1
    for i = 1:NG
        if i in ctrl.state_indices
            δ = setindex(δ, δc[iδu], i)
            iδu += 1
        else
            δ = setindex(δ, δp[iδp], i)
            iδp += 1
        end
    end

    # rigid body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # no point loads are dependent on state rates

    # save aerodynamic inputs and distributed loads
    for i = 1:NA
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_flap = flap.models[i]
        section_ctrl = LiftingLineSectionControl()
        section_models = (section_aero, section_stru, section_flap, section_ctrl)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_stru)
        Nysi = number_of_inputs(section_stru)
        Npsi = number_of_parameters(section_stru)
        Nufi = number_of_states(section_flap)
        Nyfi = number_of_inputs(section_flap)
        Npfi = number_of_parameters(section_flap)
        Nuci = number_of_states(section_ctrl)
        Nyci = number_of_inputs(section_ctrl)
        Npci = number_of_parameters(section_ctrl)

        # local section properties
        icol = stru.icol_elem[i]
        element = assembly.elements[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum

        dF_elem = SVector(dus[icol+6], dus[icol+7], dus[icol+8]) .* stru.force_scaling # internal force
        dM_elem = SVector(dus[icol+9], dus[icol+10], dus[icol+11]) .* stru.force_scaling # internal moment
        dP_elem = SVector(dus[icol+12], dus[icol+13], dus[icol+14]) .* stru.mass_scaling # linear momentum
        dH_elem = SVector(dus[icol+15], dus[icol+16], dus[icol+17]) .* stru.mass_scaling # angular momentum

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)

        # length and location of the element
        ΔL = element.L # element length
        pe = element.x + u_elem # location of element (in body frame)

        # local to body transformation
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab

        # local structural to aerodynamic transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # local freestream linear and angular velocities
        vi = -R*(CtCab'*V + GXBeam.element_linear_velocity(element, P_elem, H_elem)) # local linear freestream velocity
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem) # local angular freestream velocity

        # local freestream linear and angular accelerations
        dvi = -R*GXBeam.element_linear_velocity(element, dP_elem, dH_elem)
        dωi = R*GXBeam.element_angular_velocity(element, dP_elem, dH_elem)

        # section state variable rates
        duai = SVector{Nuai}(duas[i]) # aerodynamic state variables
        dusi = vcat(dvi, dωi) # structural state variables
        dufi = SVector{Nufi}(dufs[i]) # control surface state variables
        duci = SMatrix{1,NG}(getindex.(flap.gains, i))* dδ # controller state variables
        dui = vcat(duai, dusi, dufi, duci)

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = SMatrix{1,NG}(getindex.(flap.gains, i)) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pfi = SVector{Npfi}(pfs[i]) # control surface parameters
        pci = SVector{0,Float64}() # controller parameters
        pi = vcat(pai, psi, pfi, pci)

        # section inputs (from state rates)
        yi = get_inputs_using_state_rates(section_models, dui, ui, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)
        yfi = view(yi, Nyai + Nysi + 1 : Nyai + Nysi + Nyfi)
        yci = view(yi, Nyai + Nysi + Nyfi + 1 : Nyai + Nysi + Nyfi + Nyci)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs
        y[iyf[iyfs[i]]] .= yfi # control surface inputs

        # section aerodynamic loads (in body frame)
        fi = CtCab*R'*SVector(yi[Nyai+1], yi[Nyai+2], yi[Nyai+3])
        mi = CtCab*R'*SVector(yi[Nyai+4], yi[Nyai+5], yi[Nyai+6])

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1 : yoff+3]] = fi
        y[iys[yoff+4 : yoff+6]] = mi

        # add distributed loads to total forces and moments
        Ftot += ΔL*fi
        Mtot += cross(pe, ΔL*fi) + ΔL*mi
    end

    # save trim model inputs
    loads = SVector(Ftot[1], Ftot[2], Ftot[3], Mtot[1], Mtot[2], Mtot[3])

    for i = 1:NC
        idx = ctrl.force_indices[i]
        y[iyc[i]] = loads[idx]
    end

    return y
end

# --- convenience methods --- #

function set_additional_parameters!(padd, model::LiftingLine, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}; rho, point_conditions,
    element_loads, u, v, w, p, q, r, delta) where {NF,NG,NC,TF,TG}

    np = length(stru.icol_point)
    ne = length(stru.icol_elem)

    padd[1] = rho

    for ip = 1:np
        padd[1+6*(ip-1)+1] = point_conditions[6*(ip-1)+1]
        padd[1+6*(ip-1)+2] = point_conditions[6*(ip-1)+2]
        padd[1+6*(ip-1)+3] = point_conditions[6*(ip-1)+3]
        padd[1+6*(ip-1)+4] = point_conditions[6*(ip-1)+4]
        padd[1+6*(ip-1)+5] = point_conditions[6*(ip-1)+5]
        padd[1+6*(ip-1)+6] = point_conditions[6*(ip-1)+6]
    end

    for ie = 1:ne
        padd[1+6*np+6*(ie-1)+1] = element_loads[6*(ie-1)+1]
        padd[1+6*np+6*(ie-1)+2] = element_loads[6*(ie-1)+2]
        padd[1+6*np+6*(ie-1)+3] = element_loads[6*(ie-1)+3]
        padd[1+6*np+6*(ie-1)+4] = element_loads[6*(ie-1)+4]
        padd[1+6*np+6*(ie-1)+5] = element_loads[6*(ie-1)+5]
        padd[1+6*np+6*(ie-1)+6] = element_loads[6*(ie-1)+6]
    end

    padd[1 + 6*np + 6*ne + 1] = u
    padd[1 + 6*np + 6*ne + 2] = v
    padd[1 + 6*np + 6*ne + 3] = w
    padd[1 + 6*np + 6*ne + 4] = p
    padd[1 + 6*np + 6*ne + 5] = q
    padd[1 + 6*np + 6*ne + 6] = r
    padd[1 + 6*np + 6*ne + 6 + 1 : 1 + 6*np + 6*ne + 6 + NG - NC] = delta

    return padd
end

function separate_additional_parameters(model::LiftingLine, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, padd) where {NF,NG,NC,TF,TG}

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
    delta = view(padd, 1 + 6*np + 6*ne + 6 + 1 : 1 + 6*np + 6*ne + 6 + NG - NC)

    return (rho = rho, point_conditions = point_conditions, element_loads = element_loads,
        u = u, v = v, w = w, p = p, q = q, r = r, delta = delta)
end
