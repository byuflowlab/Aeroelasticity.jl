"""
    couple_models(aero::LiftingLine, stru::GEBT, flap::LiftingLineFlaps)

Create a coupled model using a lifting line aerodynamic model, a
geometrically exact beam theory model, and a lifting line control surface model.
This model introduces additional parameters corresponding to the freestream
air density ``\\rho_\\infty``, followed by the external loads ``F_{x,i}, F_{y,i},
F_{z,i}, M_{x,i}, M_{y,i}, M_{z,i}`` or displacements ``u_{x,i}, u_{y,i}, u_{z,i},
\\theta_{x,i}, \\theta_{y,i}, \\theta_{z,i}`` for each node, followed by the constant
distributed loads ``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i}, m_{y,i}, m_{z,i}``
applied on each beam element (excluding aerodynamic loads), followed by the body
frame linear and angular velocities ``u, v, w, p, q, r`` and the control surface
deflections ``\\delta_1, \\delta_2, \\dots, \\delta_N``.

**NOTE: When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface**
"""
function couple_models(aero::LiftingLine, stru::GEBT, flap::LiftingLineFlaps)
    return (aero, stru, flap)
end

# --- traits --- #

function number_of_additional_parameters(aero::LiftingLine, stru::GEBT, flap::LiftingLineFlaps)

    return 1 + 6*length(stru.icol_point) + 6*length(stru.icol_elem) + 6 + length(flap.gains)
end

function coupling_inplaceness(::Type{<:LiftingLine}, ::Type{<:GEBT}, ::Type{<:LiftingLineFlaps})
    return InPlace()
end

function coupling_rate_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(islinear.(coupling_rate_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function coupling_state_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(islinear.(coupling_state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

coupling_parameter_jacobian_type(::Type{<:LiftingLine}, ::Type{<:GEBT},
    ::Type{LiftingLineFlaps}) = Nonlinear()

function coupling_time_gradient_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(isempty.(coupling_time_gradient_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Empty()
    elseif all(iszero.(coupling_time_gradient_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Zeros()
    elseif all(isconstant.(coupling_time_gradient_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Constant()
    elseif all(islinear.(coupling_time_gradient_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

# --- methods --- #

function get_coupling_inputs!(y, aero::LiftingLine{NA,TA}, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, dx, x, p, t) where {NA,NF,NG,TA,TF,TG}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nuf = number_of_states(flap) # number of control surface states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyf = number_of_inputs(flap) # number of control surface inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npf = number_of_parameters(flap) # number of control surface parameters
    npadd = number_of_additional_parameters(aero, stru, flap) # number of additional parameters

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
    ipadd = npa + nps + npf + 1 : npa + nps + npf + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic parameters for each section
    iufs = state_indices(flap.models) # indices of control surface states for each section
    iyfs = input_indices(flap.models) # indices of control surface inputs for each section
    ipfs = parameter_indices(flap.models) # indices of control surface parameters for each section

    # separate state variables, inputs, and parameters
    dua = view(dx, iua) # aerodynamic rates
    ua = view(x, iua) # aerodynamic states
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dus = view(dx, ius) # structural rates
    us = view(x, ius) # structural states
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    duf = view(dx, iuf) # control surface rates
    uf = view(x, iuf) # control surface states
    yf = view(y, iyf) # control surface inputs
    pf = view(p, ipf) # control surface parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic states for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    dufs = view.(Ref(duf), iufs) # control surface rates for each section
    ufs = view.(Ref(uf), iufs) # control surface states for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

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
    δ = padd[SVector{NG}(1 + 6*npoint + 6*nelem + 6 + 1 : 1 + 6*npoint + 6*nelem + 6 + NG)]
    dδ = zero(δ)

    # body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # construct assembly from structural parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

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

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

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

        # section inputs
        yi = get_coupling_inputs(section_models, dui, ui, pi, t)

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

function set_additional_parameters!(padd, model::LiftingLine, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}; rho, point_conditions,
    element_loads, u, v, w, p, q, r, delta) where {NF,NG,TF,TG}

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
    padd[1 + 6*np + 6*ne + 6 + 1 : 1 + 6*np + 6*ne + 6 + NG] = delta

    return padd
end

function separate_additional_parameters(model::LiftingLine, stru::GEBT,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, padd) where {NF,NG,TF,TG}

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
    delta = view(padd, 1 + 6*np + 6*ne + 6 + 1 : 1 + 6*np + 6*ne + 6 + NG)

    return (rho = rho, point_conditions = point_conditions, element_loads = element_loads,
        u = u, v = v, w = w, p = p, q = q, r = r, delta = delta)
end
