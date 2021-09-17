"""
    couple_models(aero::LiftingLine, dyn::RigidBody)

Create a coupled model using a lifting line aerodynamic model and a rigid body
dynamics model.  This model introduces additional parameters corresponding to
the length, position, and orientation of each lifting line element ``L, p_e,
e_1, e_2, e_3``, followed by the freestream air density ``rho_\\infty``, rigid
body inertial properties ``m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz``, and the additional
forces/moments applied on the body ``F_x, F_y, F_z, M_x, M_y, M_z``.

**NOTE: When using this model, the local frame for each lifting line element should be
oriented with the x-axis in the chordwise direction, the y-axis in the spanwise
direction (out the right wing), and the z-axis in the airfoil normal direction**
"""
couple_models(aero::LiftingLine, dyn::RigidBody) = (aero, dyn)

# --- Traits --- #

function number_of_additional_parameters(::Type{<:LiftingLine{NA,TA}}, ::Type{<:RigidBody}) where {NA,TA}
    return 13*NA + 14
end

coupling_inplaceness(::Type{<:LiftingLine}, ::Type{<:RigidBody}) = InPlace()

function coupling_rate_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody}) where {NA,TA}
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

function coupling_state_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody}) where {NA,TA}
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

coupling_parameter_jacobian_type(::Type{<:LiftingLine}, ::Type{<:RigidBody}) = Nonlinear()

function coupling_time_gradient_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody}) where {NA,TA}
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

# --- Methods --- #

function get_coupling_inputs!(y, aero::LiftingLine{NA,TA}, dyn::RigidBody, dx, x, p, t) where {NA,TA}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nud = number_of_states(dyn) # number of rigid body states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    npadd = number_of_additional_parameters(aero, dyn) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iud = nua + 1 : nua + nud # indices of rigid body states
    iyd = nya + 1 : nya + nyd # indices of rigid body inputs
    ipd = npa + 1 : npa + npd # indices of rigid body parameters
    ipadd = npa + npd + 1 : npa + npd + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate rates, states, inputs, and parameters
    dua = view(dx, iua) # aerodynamic rates
    ua = view(x, iua) # aerodynamic states
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dud = view(dx, iud) # rigid body rates
    ud = view(x, iud) # rigid body states
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    padd = view(p, ipadd) # additional parameters
    duas = view.(Ref(dua), iuas) # aerodynamic rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic states for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body rates and states
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dud
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = ud

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

    # body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # body linear and angular acceleration
    dV = SVector(dur, dvr, dwr)
    dΩ = SVector(dpr, dqr, drr)

    # aerodynamic forces and moments
    Fa = @SVector zeros(3)
    Ma = @SVector zeros(3)
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
        offset = 13*(i-1)
        ΔL = padd[offset+1] # element length
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
        duai = SVector{Nuai}(duas[i]) # aerodynamic state variables
        dusi = vcat(dvi, dωi) # rigid body state variables
        dui = vcat(duai, dusi)

        # section states
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # rigid body state variables
        ui = vcat(uai, usi)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # rigid body parameters
        pi = vcat(pai, psi)

        # section inputs
        yi = get_coupling_inputs(section_models, dui, ui, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs

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

    return y
end

# --- Performance Overloads --- #

# TODO

# --- convenience methods --- #

function set_additional_parameters!(padd, model::LiftingLine{NA,TA},
    dyn::RigidBody; section_parameters, rho, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz,
    Fbx, Fby, Fbz, Mbx, Mby, Mbz) where {NA,TA}

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


    return padd
end

function separate_additional_parameters(model::LiftingLine{NA,TA},
    dyn::RigidBody, padd) where {NA,TA}

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
        )
end
