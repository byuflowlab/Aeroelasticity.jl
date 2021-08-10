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

function inplaceness(::Type{<:LiftingLine}, ::Type{<:RigidBody},
    ::Type{<:LiftingLineFlaps}, ::Type{<:Trim})

    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}, ::Type{<:Trim}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(isempty.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NF,NG,TF,TG}}, ::Type{<:Trim}) where {NA,NF,NG,TA,TF,TG}

    aero_model_types = TA.parameters
    flap_model_types = TF.parameters
    if all(isempty.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection),
        flap_model_types, Ref(LiftingLineSectionControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function number_of_parameters(::LiftingLine{NA,TA}, ::RigidBody,
    ::LiftingLineFlaps{NF,NG,TF,TG}, ::Trim{NC}) where {NA,NF,NG,NC,TA,TF,TG}

    return 13*NA + 14 + NG - NC
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{NA,TA}, dyn::RigidBody,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, u, p, t) where {NA,NF,NG,NC,TA,TF,TG}

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
    npadd = number_of_parameters(aero, dyn, flap, ctrl) # number of additional parameters

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

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    ud = view(u, iud) # rigid body state variables
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    uf = view(u, iuf) # rigid body state variables
    yf = view(y, iyf) # rigid body inputs
    pf = view(p, ipf) # rigid body parameters
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

    # rigid body states
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

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # rigid body state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = SMatrix{1,NG}(getindex.(flap.gains, i)) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # rigid body parameters
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

    # save trim model inputs
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

function get_input_mass_matrix!(My, aero::LiftingLine{NA,TA}, dyn::RigidBody,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, u, p, t) where {NA,NF,NG,NC,TA,TF,TG}

    # zero out mass matrix
    My .= 0

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
    npadd = number_of_parameters(aero, dyn, flap, ctrl) # number of additional parameters

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

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    ud = view(u, iud) # rigid body state variables
    pd = view(p, ipd) # rigid body parameters
    uf = view(u, iuf) # rigid body state variables
    pf = view(p, ipf) # rigid body parameters
    uc = view(u, iuc) # structural state variables
    pc = view(p, ipc) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # rigid body states
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

    # aerodynamic forces and moments
    Fa_dV = @SMatrix zeros(3, 3)
    Fa_dΩ = @SMatrix zeros(3, 3)
    Ma_dV = @SMatrix zeros(3, 3)
    Ma_dΩ = @SMatrix zeros(3, 3)
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
        ΔL = padd[offset+1] # element length
        pe = SVector(padd[offset+2], padd[offset+3], padd[offset+4]) # local frame position
        e1 = SVector(padd[offset+5], padd[offset+6], padd[offset+7]) # local frame x-axis
        e2 = SVector(padd[offset+8], padd[offset+9], padd[offset+10]) # local frame y-axis
        e3 = SVector(padd[offset+11], padd[offset+12], padd[offset+13]) # local frame z-axis
        Cab = [e1 e2 e3] # local to body frame transformation matrix
        Cba = Cab' # body to local frame transformation matrix
        vi = -Cba*(V + cross(Ω, pe)) # local linear freestream velocity
        ωi = Cba*Ω # local angular freestream velocity
        dvi_dV = -Cba
        dvi_dΩ = Cba*GXBeam.tilde(pe)
        dωi_dΩ = Cba

        # local section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # rigid body state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = SMatrix{1,NG}(getindex.(flap.gains, i)) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # local section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # rigid body parameters
        pfi = SVector{Npfi}(pfs[i]) # control surface parameters
        pci = SVector{0,Float64}() # controller parameters
        pi = vcat(pai, psi, pfi, pci)

        # local section input mass matrix
        Myi = get_input_mass_matrix(section_models, ui, pi, t)

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
        yai_dV = yai_dvi * dvi_dV
        yai_dΩ = yai_dvi * dvi_dΩ + yai_dωi * dωi_dΩ
        f_dV = f_dvi * dvi_dV
        f_dΩ = f_dvi * dvi_dΩ + f_dωi * dωi_dΩ
        m_dV = m_dvi * dvi_dV
        m_dΩ = m_dvi * dvi_dΩ + m_dωi * dωi_dΩ
        yfi_dV = yfi_dvi * dvi_dV
        yfi_dΩ = yfi_dvi * dvi_dΩ + yfi_dωi * dωi_dΩ

        # save local inputs

        # aerodynamic inputs
        My[iya[iyas[i]], iua[iuas[i]]] .= yai_duai
        My[iya[iyas[i]], iud[7:9]] .= yai_dV
        My[iya[iyas[i]], iud[10:12]] .= yai_dΩ
        My[iya[iyas[i]], iuf[iufs[i]]] .= yai_dufi

        # control surface inputs
        My[iyf[iyfs[i]], iua[iuas[i]]] .= yfi_duai
        My[iyf[iyfs[i]], iud[7:9]] .= yfi_dV
        My[iyf[iyfs[i]], iud[10:12]] .= yfi_dΩ
        My[iyf[iyfs[i]], iuf[iufs[i]]] .= yfi_dufi

        # add to global inputs

        # section loads (in body frame)
        Fi_duai = Cab*ΔL*f_duai
        Fi_dV = Cab*ΔL*f_dV
        Fi_dΩ = Cab*ΔL*f_dΩ
        Mi_duai = Cab*ΔL*m_duai
        Mi_dV = Cab*ΔL*m_dV
        Mi_dΩ = Cab*ΔL*m_dΩ

        # add loads to total aerodynamic loads
        Fa_duai = Fi_duai
        Fa_dV += Fi_dV
        Fa_dΩ += Fi_dΩ
        Ma_duai = GXBeam.tilde(pe)*Fi_duai + Mi_duai
        Ma_dV += GXBeam.tilde(pe)*Fi_dV + Mi_dV
        Ma_dΩ += GXBeam.tilde(pe)*Fi_dΩ + Mi_dΩ

        # save load mass matrix entries for this element
        My[iyd[8:10], iua[iuas[i]]] = Fa_duai
        My[iyd[11:13], iua[iuas[i]]] = Ma_duai
    end

    # save load mass matrix entries
    My[iyd[8:10], iud[7:9]] = Fa_dV
    My[iyd[8:10], iud[10:12]] = Fa_dΩ
    My[iyd[11:13], iud[7:9]] = Ma_dV
    My[iyd[11:13], iud[10:12]] = Ma_dΩ

    loads_dV = vcat(Fa_dV, Ma_dV)
    loads_dΩ = vcat(Fa_dΩ, Ma_dΩ)

    for i = 1:NC
        idx = ctrl.force_indices[i]
        My[iyc[i], iud[7:9]] = loads_dV[idx,:]
        My[iyc[i], iud[10:12]] = loads_dΩ[idx,:]
    end

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{NA,TA}, dyn::RigidBody,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, ctrl::Trim{NC}, du, u, p, t) where {NA,NF,NG,NC,TA,TF,TG}

    # initialize input vector
    models = (aero, dyn, flap, ctrl)
    T = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    y = zeros(T, Ny)

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
    npadd = number_of_parameters(aero, dyn, flap, ctrl) # number of additional parameters

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

    # separate state variables, inputs, and parameters
    dua = view(du, iua) # aerodynamic rates
    ua = view(u, iua) # aerodynamic states
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dud = view(du, iud) # rigid body rates
    ud = view(u, iud) # rigid body states
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    duf = view(du, iuf) # control surface rates
    uf = view(u, iuf) # control surface states
    yf = view(y, iyf) # control surface inputs
    pf = view(p, ipf) # control surface parameters
    duc = view(du, iuc) # controller rates
    uc = view(u, iuc) # structural states
    yc = view(y, iyc) # structural inputs
    pc = view(p, ipc) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic state variables for each section
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    dufs = view.(Ref(duf), iufs) # control surface state variables for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # rigid body states and rates
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = ud
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dud

    # controller state variables
    dδc = duc
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
    T = promote_type(eltype(du), eltype(p))
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

    # rigid body linear and angular acceleration
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

        # section state rates
        duai = SVector{Nuai}(duas[i]) # aerodynamic state variables
        dusi = vcat(dvi, dωi) # rigid body state variables
        dufi = SVector{Nufi}(dufs[i]) # control surface state variables
        duci = SMatrix{1,NG}(getindex.(flap.gains, i)) * dδ # controller state variables
        dui = vcat(duai, dusi, dufi, duci)

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # rigid body state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = SMatrix{1,NG}(getindex.(flap.gains, i)) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # section parameters
        pai = SVector{Npai}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # rigid body parameters
        pfi = SVector{Npfi}(pfs[i]) # control surface parameters
        pci = SVector{0,Float64}() # controller parameters
        pi = vcat(pai, psi, pfi, pci)

        # section inputs
        yi = get_inputs_from_state_rates(section_models, dui, ui, pi, t)

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
        Fi = Cab*ΔL*SVector(yi[Nyai+1], yi[Nyai+2], yi[Nyai+3])
        Mi = Cab*ΔL*SVector(yi[Nyai+4], yi[Nyai+5], yi[Nyai+6])

        # add loads to total aerodynamic loads
        Fa += Fi
        Ma += cross(pe, Fi) + Mi
    end

    # save rigid body inputs
    y[iyd[8]] = Fa[1]
    y[iyd[9]] = Fa[2]
    y[iyd[10]] = Fa[3]
    y[iyd[11]] = Ma[1]
    y[iyd[12]] = Ma[2]
    y[iyd[13]] = Ma[3]

    # save trim model inputs
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
