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

# --- traits --- #

inplaceness(::Type{<:LiftingLine}, ::Type{<:RigidBody}) = InPlace()

function mass_matrix_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody}) where {NA,TA}
    model_types = TA.parameters
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

function state_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:RigidBody}) where {NA,TA}
    model_types = TA.parameters
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

function number_of_parameters(::Type{<:LiftingLine{NA,TA}}, ::Type{<:RigidBody}) where {NA,TA}
    return 13*NA + 14
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{NA,TA}, dyn::RigidBody, u, p, t) where {NA,TA}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nud = number_of_states(dyn) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nyd = number_of_inputs(dyn) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    npd = number_of_parameters(dyn) # number of structural parameters
    npadd = number_of_parameters(aero, dyn) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iud = nua + 1 : nua + nud # indices of structural states
    iyd = nya + 1 : nya + nyd # indices of structural inputs
    ipd = npa + 1 : npa + npd # indices of structural parameters
    ipadd = npa + npd + 1 : npa + npd + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    ud = view(u, iud) # structural state variables
    yd = view(y, iyd) # structural inputs
    pd = view(p, ipd) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states
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

function get_input_mass_matrix!(My, aero::LiftingLine{NA,TA}, dyn::RigidBody,
    u, p, t) where {NA,TA}

    # zero out mass matrix
    My .= 0

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nud = number_of_states(dyn) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nyd = number_of_inputs(dyn) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    npd = number_of_parameters(dyn) # number of structural parameters
    npadd = number_of_parameters(aero, dyn) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iud = nua + 1 : nua + nud # indices of structural states
    iyd = nya + 1 : nya + nyd # indices of structural inputs
    ipd = npa + 1 : npa + npd # indices of structural parameters
    ipadd = npa + npd + 1 : npa + npd + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables and parameters
    ua = view(u, iua) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    ud = view(u, iud) # structural state variables
    pd = view(p, ipd) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states
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

    # aerodynamic forces and moments
    Fa_dV = @SMatrix zeros(3, 3)
    Fa_dΩ = @SMatrix zeros(3, 3)
    Ma_dV = @SMatrix zeros(3, 3)
    Ma_dΩ = @SMatrix zeros(3, 3)
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
        dvi_dV = -Cba
        dvi_dΩ = Cba*GXBeam.tilde(pe)
        dωi_dΩ = Cba

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

        # propgate derivatives using chain rule
        yai_dV = yai_dvi * dvi_dV
        yai_dΩ = yai_dvi * dvi_dΩ + yai_dωi * dωi_dΩ
        f_dV = f_dvi * dvi_dV
        f_dΩ = f_dvi * dvi_dΩ + f_dωi * dωi_dΩ
        m_dV = m_dvi * dvi_dV
        m_dΩ = m_dvi * dvi_dΩ + m_dωi * dωi_dΩ

        # save local inputs
        My[iya[iyas[i]], iua[iuas[i]]] = yai_duai
        My[iya[iyas[i]], iud[7:9]] = yai_dV
        My[iya[iyas[i]], iud[10:12]] = yai_dΩ

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

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{NA,TA}, dyn::RigidBody,
    du, u, p, t) where {NA,TA}

    # initialize input vector
    models = (aero, dyn)
    TF = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    y = zeros(TF, Ny)

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nud = number_of_states(dyn) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nyd = number_of_inputs(dyn) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    npd = number_of_parameters(dyn) # number of structural parameters
    npadd = number_of_parameters(aero, dyn) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iud = nua + 1 : nua + nud # indices of structural states
    iyd = nya + 1 : nya + nyd # indices of structural inputs
    ipd = npa + 1 : npa + npd # indices of structural parameters
    ipadd = npa + npd + 1 : npa + npd + npadd # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state rates, states, inputs, and parameters
    dua = view(du, iua) # aerodynamic states rates
    ua = view(u, iua) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dud = view(du, iud) # structural state rates
    ud = view(u, iud) # structural state variables
    yd = view(y, iyd) # structural inputs
    pd = view(p, ipd) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic state rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states and rates
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = ud
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dud

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

        # section state rates
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

        # section inputs from state rates
        yi = get_inputs_from_state_rates(section_models, dui, ui, pi, t)

        # separate inputs
        yai = view(yi, 1:Nyai)
        ysi = view(yi, Nyai + 1 : Nyai + Nysi)

        # save local inputs
        y[iya[iyas[i]]] .= yai # aerodynamic inputs

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

    return y
end
