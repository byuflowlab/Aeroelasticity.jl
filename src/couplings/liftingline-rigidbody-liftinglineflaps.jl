"""
    couple_models(aero::LiftingLine, stru::RigidBody, flap::LiftingLineFlaps)

Create an aerostructural model using a lifting line aerodynamic model, a rigid
body model, and a lifting line control surface model.  This model introduces
additional parameters corresponding to the freestream velocity components
``V_x, V_y, V_z``, air density ``\\rho``, gravitational constant ``g``,
inertial properties ``m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz``, and constant applied
forces/moments (in the body frame) ``Fx, Fy, Fz, Mx, My, Mz``, followed by the
length, position, and orientation of each lifting line element ``L, p_e, e_1,
e_2, e_3``, and the control surface deflections ``\\delta_1, \\delta_2, \\dots,
\\delta_N``.

**NOTE: When using this model, the local frame for each lifting line element should be
oriented with the x-axis in the chordwise direction, the y-axis in the spanwise
direction (out the right wing), and the z-axis in the airfoil normal direction**
"""
function couple_models(aero::LiftingLine, stru::RigidBody, flap::LiftingLineFlaps)
    return (aero, stru, flap)
end

# --- traits --- #

function inplaceness(::Type{LiftingLine}, ::Type{<:RigidBody}, ::Type{<:LiftingLineFlaps})
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{NS,TS}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NS,NF,TF}}) where {NS,NF,TS,TF}

    aero_model_types = TS.parameters
    flap_model_types = TF.parameters
    if all(isempty.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(isempty.(mass_matrix_type.(flap_model_types, Ref(LiftingLineControl))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(iszero.(mass_matrix_type.(flap_model_types, Ref(LiftingLineControl))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(isidentity.(mass_matrix_type.(flap_model_types, Ref(LiftingLineControl))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(isconstant.(mass_matrix_type.(flap_model_types, Ref(LiftingLineControl))))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types, Ref(LiftingLineSection)))) &&
        all(islinear.(mass_matrix_type.(flap_model_types, Ref(LiftingLineControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{NS,TS}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NS,NF,TF}}) where {NS,NF,TS,TF}

    aero_model_types = TS.parameters
    flap_model_types = TF.parameters
    if all(isempty.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(isempty.(state_jacobian_type.(flap_model_types, Ref(LiftingLineControl))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(iszero.(state_jacobian_type.(flap_model_types, Ref(LiftingLineControl))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(isidentity.(state_jacobian_type.(flap_model_types, Ref(LiftingLineControl))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(isconstant.(state_jacobian_type.(flap_model_types, Ref(LiftingLineControl))))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(aero_model_types, Ref(LiftingLineSection)))) &&
        all(islinear.(state_jacobian_type.(flap_model_types, Ref(LiftingLineControl))))
        return Linear()
    else
        return Nonlinear()
    end
end

function number_of_parameters(::Type{LiftingLine{NS,TS}}, ::Type{<:RigidBody},
    ::Type{LiftingLineFlaps{NS,NF,TF}}) where {NS,NF,TS,TF}

    return 18 + 13*NS + NF
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{NS,TS}, stru::RigidBody,
    flap::LiftingLineFlaps{NS,NF,TF}, u, p, t) where {NS,NF,TS,TF}

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
    npadd = number_of_parameters(aero, stru, flap) # number of additional parameters

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
    ipadd = npa + nps + npf : npa + nps + npf + npadd # indices of additional parameters
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
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # extract global state variables
    xr = us[1] # x-position
    yr = us[2] # y-position
    zr = us[3] # z-position
    ϕr = us[4] # roll
    θr = us[5] # pitch
    ψr = us[6] # yaw
    ur = us[7] # x-velocity
    vr = us[8] # y-velocity
    wr = us[9] # z-velocity
    pr = us[10] # roll rate
    qr = us[11] # pitch rate
    rr = us[12] # yaw rate

    # extract global parameters
    Vx = padd[1] # freestream x-velocity
    Vy = padd[2] # freestream y-velocity
    Vz = padd[3] # freestream z-velocity
    ρ = padd[4] # freestream air density
    g = padd[5] # gravitational constant
    m = padd[6] # rigid body mass
    Ixx = padd[7] # rigid body moment of inertia
    Iyy = padd[8] # rigid body moment of inertia
    Izz = padd[9] # rigid body moment of inertia
    Ixz = padd[10] # rigid body moment of inertia
    Ixy = padd[11] # rigid body moment of inertia
    Iyz = padd[12] # rigid body moment of inertia
    Fx = padd[13] # x-direction forces on the rigid body
    Fy = padd[14] # y-direction forces on the rigid body
    Fz = padd[15] # z-direction forces on the rigid body
    Mx = padd[16] # x-direction moment on the rigid body
    My = padd[17] # y-direction moment on the rigid body
    Mz = padd[18] # z-direction moment on the rigid body
    δ = padd[SVector{NF}(18 + 13*NS + 1 : 18 + 13*NS + NF)] # control surface deflections

    # rigid body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - V

    # section forces and moments
    Fa = @SVector zeros(3)
    Ma = @SVector zeros(3)
    for i = 1:N
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_flap = flap.models[i]
        section_ctrl = LiftingLineControl()
        section_models = (section_aero, section_stru, section_flap, section_ctrl)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_aero)
        Nysi = number_of_inputs(section_aero)
        Npsi = number_of_parameters(section_aero)
        Nufi = number_of_states(section_flap)
        Nyfi = number_of_inputs(section_flap)
        Npfi = number_of_parameters(section_flap)
        Nuci = number_of_states(section_aero)
        Nyci = number_of_inputs(section_aero)
        Npci = number_of_parameters(section_aero)

        # local section properties
        offset = 18 + 13*(i-1)
        ΔL = pc[offset+1] # section spanwise length
        pe = SVector(pc[offset+2], pc[offset+3], pc[offset+4]) # local frame position
        e1 = SVector(pc[offset+5], pc[offset+6], pc[offset+7]) # local frame x-axis
        e2 = SVector(pc[offset+8], pc[offset+9], pc[offset+10]) # local frame y-axis
        e3 = SVector(pc[offset+11], pc[offset+12], pc[offset+13]) # local frame z-axis
        Cab = [e1 e2 e3] # local to body frame transformation matrix
        Cba = Cab' # body to local frame transformation matrix
        vi = Cba*Vinf - Cba*cross(Ω, pe) # local linear freestream velocity
        ωi = Cba*Ω # local angular freestream velocity

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = getindex.(flap, i) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # section parameters
        pai = SVector{Nupi}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pfi = SVector{Nufi}(pfs[i]) # control surface parameters
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
        y[ifa[iyfs[i]]] .= yfi # control surface inputs

        # add to global inputs

        # section aerodynamic loads (in body frame)
        Fi = Cab*ΔL*SVector(yi[Nyai+1], yi[Nyai+2], yi[Nyai+3])
        Mi = Cab*ΔL*SVector(yi[Nyai+4], yi[Nyai+5], yi[Nyai+6])

        # add loads to total aerodynamic loads
        Fa += Fi
        Ma += cross(pe, Fi) + Mi
    end

    # gravitational forces and moments
    sϕ, cϕ = sincos(ϕr)
    sθ, cθ = sincos(θr)
    Fg = m*g*SVector(-sθ, cθ*sϕ, cθ*cϕ)
    Mg = @SVector zeros(3)  # no gravitational moment

    # save rigid body inputs
    y[iys[1]] = m
    y[iys[2]] = Ixx
    y[iys[3]] = Iyy
    y[iys[4]] = Izz
    y[iys[5]] = Ixz
    y[iys[6]] = Ixy
    y[iys[7]] = Iyz
    y[iys[8]] = Fa[1] + Fg[1] + Fx
    y[iys[9]] = Fa[2] + Fg[2] + Fy
    y[iys[10]] = Fa[3] + Fg[3] + Fz
    y[iys[11]] = Ma[1] + Mg[1] + Mx
    y[iys[12]] = Ma[2] + Mg[2] + My
    y[iys[13]] = Ma[3] + Mg[3] + Mz

    return y
end

function get_input_mass_matrix!(M, aero::LiftingLine{NS,TS}, stru::RigidBody,
    flap::LiftingLineFlaps{NS,NF,TF}, u, p, t) where {NS,NF,TS,TF}

    # zero out mass matrix
    M .= 0

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
    npadd = number_of_parameters(aero, stru, flap) # number of additional parameters

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
    ipadd = npa + nps + npf : npa + nps + npf + npadd # indices of additional parameters
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
    padd = view(p, ipadd) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # extract global state variables
    xr = us[1] # x-position
    yr = us[2] # y-position
    zr = us[3] # z-position
    ϕr = us[4] # roll
    θr = us[5] # pitch
    ψr = us[6] # yaw
    ur = us[7] # x-velocity
    vr = us[8] # y-velocity
    wr = us[9] # z-velocity
    pr = us[10] # roll rate
    qr = us[11] # pitch rate
    rr = us[12] # yaw rate

    # extract global parameters
    Vx = padd[1] # freestream x-velocity
    Vy = padd[2] # freestream y-velocity
    Vz = padd[3] # freestream z-velocity
    ρ = padd[4] # freestream air density
    g = padd[5] # gravitational constant
    m = padd[6] # rigid body mass
    Ixx = padd[7] # rigid body moment of inertia
    Iyy = padd[8] # rigid body moment of inertia
    Izz = padd[9] # rigid body moment of inertia
    Ixz = padd[10] # rigid body moment of inertia
    Ixy = padd[11] # rigid body moment of inertia
    Iyz = padd[12] # rigid body moment of inertia
    Fx = padd[13] # x-direction forces on the rigid body
    Fy = padd[14] # y-direction forces on the rigid body
    Fz = padd[15] # z-direction forces on the rigid body
    Mx = padd[16] # x-direction moment on the rigid body
    My = padd[17] # y-direction moment on the rigid body
    Mz = padd[18] # z-direction moment on the rigid body
    δ = padd[SVector{NF}(18 + 13*NS + 1 : 18 + 13*NS + NF)] # control surface deflections

    # rigid body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - V

    # aerodynamic forces and moments
    Fa_dV = @SMatrix zeros(3, 3)
    Fa_dΩ = @SMatrix zeros(3, 3)
    Ma_dV = @SMatrix zeros(3, 3)
    Ma_dΩ = @SMatrix zeros(3, 3)
    for i = 1:N
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_flap = flap.models[i]
        section_ctrl = LiftingLineControl()
        section_models = (section_aero, section_stru, section_flap, section_ctrl)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_aero)
        Nysi = number_of_inputs(section_aero)
        Npsi = number_of_parameters(section_aero)
        Nufi = number_of_states(section_flap)
        Nyfi = number_of_inputs(section_flap)
        Npfi = number_of_parameters(section_flap)
        Nuci = number_of_states(section_aero)
        Nyci = number_of_inputs(section_aero)
        Npci = number_of_parameters(section_aero)

        # local section properties
        offset = npa + nps + 18 + 13*(i-1)
        ΔL = p[offset+1] # element length
        pe = SVector(p[offset+2], p[offset+3], p[offset+4]) # local frame position
        e1 = SVector(p[offset+5], p[offset+6], p[offset+7]) # local frame x-axis
        e2 = SVector(p[offset+8], p[offset+9], p[offset+10]) # local frame y-axis
        e3 = SVector(p[offset+11], p[offset+12], p[offset+13]) # local frame z-axis
        Cab = [e1 e2 e3] # local to body frame transformation matrix
        Cba = Cab' # body to local frame transformation matrix
        vi = Cba*Vinf - Cba*cross(Ω, pe) # local linear freestream velocity
        ωi = Cba*Ω # local angular freestream velocity
        dvi_dV = -Cba
        dvi_dΩ = Cba*GXBeam.tilde(pe)
        dωi_dΩ = Cba

        # local section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = getindex.(flap, i) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # local section parameters
        pai = SVector{Nupi}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pfi = SVector{Nufi}(pfs[i]) # control surface parameters
        pci = SVector{0,Float64}() # controller parameters
        pi = vcat(pai, psi, pfi, pci)

        # local section input mass matrix
        Mi = get_input_mass_matrix(model, LiftingLineSection(), ui, pi, t)

        # separate into component mass matrices
        yai_duai = SMatrix{Nyai,Nuai}(view(Mi, 1:Nyai, 1:Nuai))
        yai_dvi = SMatrix{Nyai,3}(view(Mi, 1:Nyai, Nuai+1:Nuai+3))
        yai_dωi = SMatrix{Nyai,3}(view(Mi, 1:Nyai, Nuai+4:Nuai+6))
        yai_dufi = SMatrix{Nyai,Nufi}(view(Mi, 1:Nyai, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        yai_duci = SMatrix{Nyai,Nuci}(view(Mi, 1:Nyai, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        f_duai = SMatrix{3,Nuai}(view(Mi, Nyai+1:Nyai+3, 1:Nuai))
        f_dvi = SMatrix{3,3}(view(Mi, Nyai+1:Nyai+3, Nuai+1:Nuai+3))
        f_dωi = SMatrix{3,3}(view(Mi, Nyai+1:Nyai+3, Nuai+4:Nuai+6))
        f_dufi = SMatrix{3,Nufi}(view(Mi, Nyai+1:Nyai+3, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        f_duci = SMatrix{3,Nuci}(view(Mi, Nyai+1:Nyai+3, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        m_duai = SMatrix{3,Nuai}(view(Mi, Nyai+4:Nyai+6, 1:Nuai))
        m_dvi = SMatrix{3,3}(view(Mi, Nyai+4:Nyai+6, Nuai+1:Nuai+3))
        m_dωi = SMatrix{3,3}(view(Mi, Nyai+4:Nyai+6, Nuai+4:Nuai+6))
        m_dufi = SMatrix{3,Nufi}(view(Mi, Nyai+4:Nyai+6, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        m_duci = SMatrix{3,Nuci}(view(Mi, Nyai+4:Nyai+6, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        yfi_duai = SMatrix{Nyfi,Nuai}(view(Mi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, 1:Nuai))
        yfi_dvi = SMatrix{Nyfi,3}(view(Mi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+1:Nuai+3))
        yfi_dωi = SMatrix{Nyfi,3}(view(Mi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+4:Nuai+6))
        yfi_dufi = SMatrix{Nyfi,Nufi}(view(Mi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        yfi_duci = SMatrix{Nyfi,Nuci}(view(Mi, Nyai+Nysi+1:Nyai+Nysi+Nyfi, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))
        yci_duai = SMatrix{Nyci,Nuai}(view(Mi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, 1:Nuai))
        yci_dvi = SMatrix{Nyci,3}(view(Mi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+1:Nuai+3))
        yci_dωi = SMatrix{Nyci,3}(view(Mi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+4:Nuai+6))
        yci_dufi = SMatrix{Nyci,Nufi}(view(Mi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+Nusi+1:Nuai+Nusi+Nufi))
        yci_duci = SMatrix{Nyci,Nuci}(view(Mi, Nyai+Nysi+Nyfi+1:Nyai+Nysi+Nyfi+Nyci, Nuai+Nusi+Nufi+1:Nuai+Nusi+Nufi+Nuci))

        # propagate derivatives using chain rule
        yai_dV = yai_dvi * dvi_dV
        yai_dΩ = yai_dvi * dvi_dΩ + yai_dωi * dωi_dΩ
        f_dV = f_dvi * dvi_dV
        f_dΩ = f_dvi * dvi_dΩ + f_dωi * dωi_dΩ
        m_dV = m_dvi * dvi_dV
        m_dΩ = m_dvi * dvi_dΩ + m_dωi * dωi_dΩ
        yfi_dV = yfi_dvi * dvi_dV
        yfi_dΩ = yfi_dvi * dvi_dΩ + yai_dωi * dωi_dΩ

        # save local inputs

        # aerodynamic inputs
        M[iya[iyas[i]], iua[iuas[i]]] .= yai_duai
        M[iya[iyas[i]], ius[7:9]] .= yai_dV
        M[iya[iyas[i]], ius[10:12]] .= yai_dΩ
        M[iya[iyas[i]], iuf[iufs[i]]] .= yai_dufi

        # control surface inputs
        M[ifa[iyfs[i]], iua[iuas[i]]] .= yfi_duai
        M[ifa[iyfs[i]], ius[7:9]] .= yfi_dV
        M[ifa[iyfs[i]], ius[10:12]] .= yfi_dΩ
        M[ifa[iyfs[i]], iuf[iufs[i]]] .= yfi_dufi

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
        M[iys[8:10], iua[iuas[i]]] = Fa_duai
        M[iys[11:13], iua[iuas[i]]] = Ma_duai
    end
    # save load mass matrix entries
    M[iys[8:10], ius[7:9]] = Fa_dV
    M[iys[8:10], ius[10:12]] = Fa_dΩ
    M[iys[11:13], ius[7:9]] = Ma_dV
    M[iys[11:13], ius[10:12]] = Ma_dΩ

    return M
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{NS,TS}, stru::RigidBody,
    flap::LiftingLineFlaps{NS,NF,TF}, du, u, p, t) where {NS,NF,TS,TF}

    # initialize input vector
    models = (aero, stru, flap)
    T = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    y = zeros(T, Ny)

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
    npadd = number_of_parameters(aero, stru, flap) # number of additional parameters

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
    ipadd = npa + nps + npf : npa + nps + npf + npadd # indices of additional parameters
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
    padd = view(p, ipadd) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic state variables for each section
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    dufs = view.(Ref(duf), iufs) # control surface state variables for each section
    ufs = view.(Ref(uf), iufs) # control surface state variables for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # extract global state variables
    xr = us[1] # x-position
    yr = us[2] # y-position
    zr = us[3] # z-position
    ϕr = us[4] # roll
    θr = us[5] # pitch
    ψr = us[6] # yaw
    ur = us[7] # x-velocity
    vr = us[8] # y-velocity
    wr = us[9] # z-velocity
    pr = us[10] # roll rate
    qr = us[11] # pitch rate
    rr = us[12] # yaw rate

    # extract global state rates
    dxr = dus[1] # x-position
    dyr = dus[2] # y-position
    dzr = dus[3] # z-position
    dϕr = dus[4] # roll
    dθr = dus[5] # pitch
    dψr = dus[6] # yaw
    dur = dus[7] # x-velocity
    dvr = dus[8] # y-velocity
    dwr = dus[9] # z-velocity
    dpr = dus[10] # roll rate
    dqr = dus[11] # pitch rate
    drr = dus[12] # yaw rate

    # extract global parameters
    Vx = padd[1] # freestream x-velocity
    Vy = padd[2] # freestream y-velocity
    Vz = padd[3] # freestream z-velocity
    ρ = padd[4] # freestream air density
    g = padd[5] # gravitational constant
    m = padd[6] # rigid body mass
    Ixx = padd[7] # rigid body moment of inertia
    Iyy = padd[8] # rigid body moment of inertia
    Izz = padd[9] # rigid body moment of inertia
    Ixz = padd[10] # rigid body moment of inertia
    Ixy = padd[11] # rigid body moment of inertia
    Iyz = padd[12] # rigid body moment of inertia
    Fx = padd[13] # x-direction forces on the rigid body
    Fy = padd[14] # y-direction forces on the rigid body
    Fz = padd[15] # z-direction forces on the rigid body
    Mx = padd[16] # x-direction moment on the rigid body
    My = padd[17] # y-direction moment on the rigid body
    Mz = padd[18] # z-direction moment on the rigid body
    δ = padd[SVector{NF}(18 + 13*NS + 1 : 18 + 13*NS + NF)] # control surface deflections

    # rigid body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # rigid body linear and angular acceleration
    dV = SVector(dur, dvr, dwr)
    dΩ = SVector(dpr, dqr, drr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - V
    dVinf = -dV

    # section forces and moments
    Fa = @SVector zeros(3)
    Ma = @SVector zeros(3)
    for i = 1:N
        # models for this section
        section_aero = aero.models[i]
        section_stru = LiftingLineSection()
        section_flap = flap.models[i]
        section_ctrl = LiftingLineControl()
        section_models = (section_aero, section_stru, section_flap, section_ctrl)

        # model dimensions for this section
        Nuai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nusi = number_of_states(section_aero)
        Nysi = number_of_inputs(section_aero)
        Npsi = number_of_parameters(section_aero)
        Nufi = number_of_states(section_flap)
        Nyfi = number_of_inputs(section_flap)
        Npfi = number_of_parameters(section_flap)
        Nuci = number_of_states(section_aero)
        Nyci = number_of_inputs(section_aero)
        Npci = number_of_parameters(section_aero)

        # local section properties
        offset = 18 + 13*(i-1)
        ΔL = pc[offset+1] # section spanwise length
        pe = SVector(pc[offset+2], pc[offset+3], pc[offset+4]) # local frame position
        e1 = SVector(pc[offset+5], pc[offset+6], pc[offset+7]) # local frame x-axis
        e2 = SVector(pc[offset+8], pc[offset+9], pc[offset+10]) # local frame y-axis
        e3 = SVector(pc[offset+11], pc[offset+12], pc[offset+13]) # local frame z-axis
        Cab = [e1 e2 e3] # local to body frame transformation matrix
        Cba = Cab' # body to local frame transformation matrix
        vi = Cba*Vinf - Cba*cross(Ω, pe) # local linear freestream velocity
        ωi = Cba*Ω # local angular freestream velocity
        dvi = Cba*dVinf - Cba*cross(dΩ, pe) # local linear freestream acceleration
        dωi = Cba*dΩ # local angular freestream acceleration

        # section state variable rates
        duai = SVector{Nuai}(duas[i]) # aerodynamic state variables
        dusi = vcat(dvi, dωi) # structural state variables
        dufi = SVector{Nufi}(dufs[i]) # control surface state variables
        duci = getindex.(flap, i) * dδ # controller state variables
        dui = vcat(duai, dusi, dufi, duci)

        # section state variables
        uai = SVector{Nuai}(uas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ufi = SVector{Nufi}(ufs[i]) # control surface state variables
        uci = getindex.(flap, i) * δ # controller state variables
        ui = vcat(uai, usi, ufi, uci)

        # section parameters
        pai = SVector{Nupi}(pas[i]) # aerodynamic parameters
        psi = SVector(ρ) # structural parameters
        pfi = SVector{Nufi}(pfs[i]) # control surface parameters
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
        y[ifa[iyfs[i]]] .= yfi # control surface inputs

        # add to global inputs

        # section aerodynamic loads (in body frame)
        Fi = Cab*ΔL*SVector(yi[Nyai+1], yi[Nyai+2], yi[Nyai+3])
        Mi = Cab*ΔL*SVector(yi[Nyai+4], yi[Nyai+5], yi[Nyai+6])

        # add loads to total aerodynamic loads
        Fa += Fi
        Ma += cross(pe, Fi) + Mi
    end

    # save total forces and moments
    y[iys[8:10]] = Fa
    y[iys[11:13]] = Ma

    return y
end
