"""
    couple_models(aero::LiftingLine, dyn::RigidBody, flap::LiftingLineFlaps, ctrl::Trim)

Create a coupled model using a lifting line aerodynamic model, a rigid
body dynamics model, a lifting line control surface model, and a controller which
maintains trimmed operating conditions.  This model introduces additional
parameters corresponding to the length, position, and orientation of each lifting
line element ``L, p_e, e_1, e_2, e_3``, followed by the freestream air density
``rho_\\infty``, rigid body inertial properties ``m, Ixx, Iyy, Izz, Ixz, Ixy,
Iyz``, additional forces/moments applied on the body ``F_x, F_y, F_z, M_x, M_y,
M_z``, and the control inputs which are not defined by the controller.

**NOTE: When using this model, the local frame for each lifting line element should be
oriented with the x-axis in the chordwise direction, the y-axis in the spanwise
direction (out the right wing), and the z-axis in the airfoil normal direction**
"""
function couple_models(aero::LiftingLine, dyn::RigidBody, flap::LiftingLineFlaps,
    ctrl::Trim)

    return (aero, dyn, flap, ctrl)
end

# --- traits --- #

function number_of_additional_parameters(::LiftingLine{NA,TA}, ::RigidBody,
    ::LiftingLineFlaps{NF,NG,TF,TG}, ::Trim{NC}) where {NA,NF,NG,NC,TA,TF,TG}

    return 13*NA + 14 + NG - NC
end

function coupling_inplaceness(::Type{<:LiftingLine}, ::Type{<:RigidBody},
    ::Type{<:LiftingLineFlaps}, ::Type{<:Trim})

    return InPlace()
end

function coupling_rate_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}, ::Type{<:Trim}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(islinear.(coupling_rate_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function coupling_state_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}, ::Type{<:Trim}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(islinear.(coupling_state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function coupling_parameter_jacobian_type(::Type{<:LiftingLine}, ::Type{<:RigidBody},
    ::Type{<:LiftingLineFlaps}, ::Type{<:Trim})

    return Nonlinear()
end

function coupling_time_gradient_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}, ::Type{<:Trim}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(isempty.(coupling_time_gradient_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Empty()
    elseif all(iszero.(coupling_time_gradient_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Zeros()
    elseif all(isinvariant.(coupling_time_gradient_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Invariant()
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

function get_coupling_inputs!(y, aero::LiftingLine{NA,TA}, dyn::RigidBody,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, dx, x, p, t) where {NA,NF,NG,NC,TA,TF,TG}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nud = number_of_states(dyn) # number of rigid body states
    nuf = number_of_states(flap) # number of control surface states
    nuc = number_of_states(ctrl) # number of control surface states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    nyf = number_of_inputs(flap) # number of control surface inputs
    nyc = number_of_inputs(ctrl) # number of control surface inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    npf = number_of_parameters(flap) # number of control surface parameters
    npc = number_of_parameters(ctrl) # number of control surface parameters
    npadd = number_of_additional_parameters(aero, dyn, flap, ctrl) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iud = nua + 1 : nua + nud # indices of rigid body states
    iyd = nya + 1 : nya + nyd # indices of rigid body inputs
    ipd = npa + 1 : npa + npd # indices of rigid body parameters
    iuf = nua + nud + 1 : nua + nud + nuf # indices of control surface states
    iyf = nya + nyd + 1 : nya + nyd + nyf # indices of control surface inputs
    ipf = npa + npd + 1 : npa + npd + npf # indices of control surface parameters
    iuc = nua + nud + nuf + 1 : nua + nud + nuf + nuc # indices of controller states
    iyc = nya + nyd + nyf + 1 : nya + nyd + nyf + nyc # indices of controller inputs
    ipc = npa + npd + npf + 1 : npa + npd + npf + npc # indices of controller parameters
    ipadd = npa + npd + npf + npc + 1 : npa + npd + npf + npc + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic parameters for each section
    iufs = state_indices(flap.models) # indices of control surface states for each section
    iyfs = input_indices(flap.models) # indices of control surface inputs for each section
    ipfs = parameter_indices(flap.models) # indices of control surface parameters for each section

    # separate rates, states, inputs, and parameters
    dua = view(dx, iua) # aerodynamic rates
    ua = view(x, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dud = view(dx, iud) # rigid body rates
    ud = view(x, iud) # rigid body state variables
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    duf = view(dx, iuf) # control surface rates
    uf = view(x, iuf) # control surface states
    yf = view(y, iyf) # control surface inputs
    pf = view(p, ipf) # control surface parameters
    duc = view(dx, iuc) # controller rates
    uc = view(x, iuc) # controller states
    yc = view(y, iyc) # controller inputs
    pc = view(p, ipc) # controller parameters
    padd = view(p, ipadd) # additional parameters
    duas = view.(Ref(dua), iuas) # aerodynamic rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic states for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    dufs = view.(Ref(duf), iufs) # control surface rates for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # rigid body rates and states
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dud
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = ud

    # controller state variables
    δc = uc

    # global parameters
    ρ = padd[13*NA+1]
    m = padd[13*NA+2]
    Ixx = padd[13*NA+3]
    Iyy = padd[13*NA+4]
    Izz = padd[13*NA+5]
    Ixz = padd[13*NA+6]
    Ixy = padd[13*NA+7]
    Iyz = padd[13*NA+8]
    Fbx = padd[13*NA+9]
    Fby = padd[13*NA+10]
    Fbz = padd[13*NA+11]
    Mbx = padd[13*NA+12]
    Mby = padd[13*NA+13]
    Mbz = padd[13*NA+14]
    δp = padd[SVector{NG-NC}(13*NA + 14 + 1 : 13*NA + 14 + NG - NC)]

    # commanded control deflections
    T = promote_type(eltype(x), eltype(p))
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
    dδ = zero(δ)

    # body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # body linear and angular acceleration
    dV = SVector(dur, dvr, dwr)
    dΩ = SVector(dpr, dqr, drr)

    # section forces and moments
    Fa = @SVector zeros(3)
    Ma = @SVector zeros(3)
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
        offset = 13*(i-1)
        ΔL = padd[offset+1] # section spanwise length
        pe = SVector(padd[offset+2], padd[offset+3], padd[offset+4]) # local frame position
        e1 = SVector(padd[offset+5], padd[offset+6], padd[offset+7]) # local frame x-axis
        e2 = SVector(padd[offset+8], padd[offset+9], padd[offset+10]) # local frame y-axis
        e3 = SVector(padd[offset+11], padd[offset+12], padd[offset+13]) # local frame z-axis
        Cab = [e1 e2 e3] # local to body frame transformation matrix
        Cba = Cab' # body to local frame transformation matrix
        vi = -Cba*(V + cross(Ω, pe)) # local linear freestream velocity
        ωi = Cba*Ω # local angular freestream velocity
        dvi =-Cba*(dV + cross(dΩ, pe)) # local linear freestream acceleration
        dωi = Cba*dΩ # local angular freestream acceleration

        # section rates
        dxai = SVector{Nuai}(duas[i]) # aerodynamic state variables
        dxsi = vcat(dvi, dωi) # rigid body state variables
        dxfi = SVector{Nufi}(dufs[i]) # control surface state variables
        dxci = SMatrix{1,NG}(getindex.(flap.gains, i)) * dδ # controller state variables
        dxi = vcat(dxai, dxsi, dxfi, dxci)

        # section states
        xai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        xsi = vcat(vi, ωi) # rigid body state variables
        xfi = SVector{Nufi}(ufs[i]) # control surface state variables
        xci = SMatrix{1,NG}(getindex.(flap.gains, i)) * δ # controller state variables
        xi = vcat(xai, xsi, xfi, xci)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # rigid body parameters
        pfi = SVector{Npfi}(pfs[i]) # control surface parameters
        pci = SVector{0,Float64}() # controller parameters
        pi = vcat(pai, psi, pfi, pci)

        # section inputs
        yi = get_coupling_inputs(section_models, dxi, xi, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)
        yfi = view(yi, Nyai + Nysi + 1 : Nyai + Nysi + Nyfi)
        yci = view(yi, Nyai + Nysi + Nyfi + 1 : Nyai + Nysi + Nyfi + Nyci)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs
        y[iyf[iyfs[i]]] .= yfi # control surface inputs

        # add to global inputs

        # section aerodynamic loads (in body frame)
        Fi = Cab*ΔL*SVector(ysi[1], ysi[2], ysi[3])
        Mi = Cab*ΔL*SVector(ysi[4], ysi[5], ysi[6])

        # add loads to total aerodynamic loads
        Fa += Fi
        Ma += cross(pe, Fi) + Mi
    end

    # save rigid body inputs
    y[iyd[1]] = m
    y[iyd[2]] = Ixx
    y[iyd[3]] = Iyy
    y[iyd[4]] = Izz
    y[iyd[5]] = Ixz
    y[iyd[6]] = Ixy
    y[iyd[7]] = Iyz
    y[iyd[8]] = Fa[1] + Fbx
    y[iyd[9]] = Fa[2] + Fby
    y[iyd[10]] = Fa[3] + Fbz
    y[iyd[11]] = Ma[1] + Mbx
    y[iyd[12]] = Ma[2] + Mby
    y[iyd[13]] = Ma[3] + Mbz

    # save controller inputs
    loads = SVector(
        Fa[1] + Fbx,
        Fa[2] + Fby,
        Fa[3] + Fbz,
        Ma[1] + Mbx,
        Ma[2] + Mby,
        Ma[3] + Mbz)

    for i = 1:NC
        idx = ctrl.force_indices[i]
        y[iyc[i]] = loads[idx]
    end

    return y
end

# --- Performance Overloads --- #

# TODO: Rate, state, input, and parameter jacobians

# --- Convenience Methods --- #

function set_additional_parameters!(padd, model::LiftingLine{NA,TA},
    dyn::RigidBody, flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC};
    section_parameters, rho, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fbx, Fby, Fbz,
    Mbx, Mby, Mbz, delta) where {NA,NF,NG,NC,TA,TF,TG}

    for i = 1:NA
        padd[13*(i-1)+1] = section_parameters[i].L
        padd[13*(i-1)+2 : 13*(i-1)+4] = section_parameters[i].pe
        padd[13*(i-1)+5 : 13*(i-1)+7] = section_parameters[i].e1
        padd[13*(i-1)+8 : 13*(i-1)+10] = section_parameters[i].e2
        padd[13*(i-1)+11 : 13*(i-1)+13] = section_parameters[i].e3
    end

    padd[13*NA+1] = rho
    padd[13*NA+2] = m
    padd[13*NA+3] = Ixx
    padd[13*NA+4] = Iyy
    padd[13*NA+5] = Izz
    padd[13*NA+6] = Ixz
    padd[13*NA+7] = Ixy
    padd[13*NA+8] = Iyz
    padd[13*NA+9] = Fbx
    padd[13*NA+10] = Fby
    padd[13*NA+11] = Fbz
    padd[13*NA+12] = Mbx
    padd[13*NA+13] = Mby
    padd[13*NA+14] = Mbz
    padd[13*NA+14+1 : 13*NA+14+NG-NC] = delta

    return padd
end

function separate_additional_parameters(model::LiftingLine{NA,TA},
    dyn::RigidBody, flap::LiftingLineFlaps{NF,NG,TF,TG}, trim::Trim{NC},
    padd) where {NA,NF,NG,NC,TA,TF,TG,TI}

    # section parameters
    return (
        section_parameters = ntuple(i -> (
            L = padd[13*(i-1)+1], # element length
            pe = view(padd, 13*(i-1)+2 : 13*(i-1)+4), # local frame position
            e1 = view(padd, 13*(i-1)+5 : 13*(i-1)+7), # local frame x-axis
            e2 = view(padd, 13*(i-1)+8 : 13*(i-1)+10), # local frame y-axis
            e3 = view(padd, 13*(i-1)+11 : 13*(i-1)+13), # local frame z-axis
            ), NA),
        rho = padd[13*NA+1],
        m = padd[13*NA+2],
        Ixx = padd[13*NA+3],
        Iyy = padd[13*NA+4],
        Izz = padd[13*NA+5],
        Ixz = padd[13*NA+6],
        Ixy = padd[13*NA+7],
        Iyz = padd[13*NA+8],
        Fbx = padd[13*NA+9],
        Fby = padd[13*NA+10],
        Fbz = padd[13*NA+11],
        Mbx = padd[13*NA+12],
        Mby = padd[13*NA+13],
        Mbz = padd[13*NA+14],
        delta = view(padd, 13*NA+14+1 : 13*NA+14+NG-NC),
        )
end
