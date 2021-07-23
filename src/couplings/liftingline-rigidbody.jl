"""
    couple_models(aero::LiftingLine, stru::RigidBody)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a rigid body model.  This model introduces additional parameters
corresponding to the freestream velocity components ``V_x, V_y, V_z``,
air density ``\\rho``, gravitational constant ``g``, inertial properties ``m,
Ixx, Iyy, Izz, Ixz, Ixy, Iyz``, and constant applied forces/moments (in the body
frame) ``Fx, Fy, Fz, Mx, My, Mz``, followed by the length, position, and
orientation of each lifting line element ``L, p_e, e_1, e_2, e_3``.

**NOTE: When using this model, the local frame for each lifting line element should be
oriented with the x-axis in the chordwise direction, the y-axis in the spanwise
direction (out the right wing), and the z-axis in the airfoil normal direction**
"""
couple_models(aero::LiftingLine, stru::RigidBody) = (aero, stru)

# --- traits --- #

function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
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

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
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

function number_of_parameters(::Type{<:LiftingLine{N,T}}, ::Type{<:RigidBody}) where {N,T}
    return 13*N+18
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{N,T}, stru::RigidBody, u, p, t) where {N,T}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npc = number_of_parameters(aero, stru) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    ipc = npa + nps + 1 : npa + nps + npc # indices of additional parameters
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
    pc = view(p, ipc) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = us

    # global parameters (first 18 additional parameters)
    Vx, Vy, Vz, ρ, g, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz = pc

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - SVector(ur, vr, wr)

    # body angular velocity
    Ω = SVector(pr, qr, rr)

    # aerodynamic forces and moments
    Fa = @SVector zeros(3)
    Ma = @SVector zeros(3)
    for i = 1:N
        # section aerodynamic model
        model = aero.models[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)
        # calculate section properties
        offset = 18 + 13*(i-1)
        ΔL = pc[offset+1] # element length
        pe = SVector(pc[offset+2], pc[offset+3], pc[offset+4]) # local frame position
        e1 = SVector(pc[offset+5], pc[offset+6], pc[offset+7]) # local frame x-axis
        e2 = SVector(pc[offset+8], pc[offset+9], pc[offset+10]) # local frame y-axis
        e3 = SVector(pc[offset+11], pc[offset+12], pc[offset+13]) # local frame z-axis
        Cab = [e1 e2 e3] # local to body frame transformation matrix
        Cba = Cab' # body to local frame transformation matrix
        vi = Cba*Vinf - Cba*cross(Ω, pe) # local linear freestream velocity
        ωi = Cba*Ω # local angular freestream velocity
        uai = SVector{Nui}(uas[i]) # aerodynamic state variables
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters
        # calculate aerodynamic inputs and section loads
        ui = vcat(uai, vi, ωi) # section state variables
        pi = vcat(pai, ρ) # section parameters
        yi = get_inputs(model, LiftingLineSection(), ui, pi, t)
        # save aerodynamic inputs
        y[iya[iyas[i]]] .= view(yi, 1:Nyi)
        # section aerodynamic loads (in body frame)
        Fi = Cab*ΔL*SVector(yi[Nyi+1], yi[Nyi+2], yi[Nyi+3])
        Mi = Cab*ΔL*SVector(yi[Nyi+4], yi[Nyi+5], yi[Nyi+6])
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

function get_input_mass_matrix!(M, aero::LiftingLine{N,T}, stru::RigidBody,
    u, p, t) where {N,T}

    # zero out mass matrix
    M .= 0

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npc = number_of_parameters(aero, stru) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    ipc = npa + nps + 1 : npa + nps + npc # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables and parameters
    ua = view(u, iua) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    us = view(u, ius) # structural state variables
    ps = view(p, ips) # structural parameters
    pc = view(p, ipc) # additional parameters for coupled model
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = us

    # global parameters (first 18 additional parameters)
    Vx, Vy, Vz, ρ, g, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz = pc

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - SVector(ur, vr, wr)

    # body angular velocity
    Ω = SVector(pr, qr, rr)

    # aerodynamic forces and moments
    Fa_dV = @SMatrix zeros(3, 3)
    Fa_dΩ = @SMatrix zeros(3, 3)
    Ma_dV = @SMatrix zeros(3, 3)
    Ma_dΩ = @SMatrix zeros(3, 3)
    for i = 1:N
        # section aerodynamic model
        model = aero.models[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)
        # calculate section properties
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
        uai = SVector{Nui}(uas[i]) # aerodynamic state variables
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters
        dvi_dV = -Cba
        dvi_dΩ = Cba*GXBeam.tilde(pe)
        dωi_dΩ = Cba
        # calculate lifting line section mass matrices
        ui = vcat(uai, vi, ωi) # section state variables
        pi = vcat(pai, ρ) # section parameters
        Mi = get_input_mass_matrix(model, LiftingLineSection(), ui, pi, t)
        # separate into component mass matrices
        d_dλi = SMatrix{Nyi,Nui}(view(Mi, 1:Nyi, 1:Nui))
        d_dvi = SMatrix{Nyi,3}(view(Mi, 1:Nyi, Nui+1:Nui+3))
        d_dωi = SMatrix{Nyi,3}(view(Mi, 1:Nyi, Nui+4:Nui+6))
        f_dλi = view(Mi, Nyi+1:Nyi+3, 1:Nui)
        f_dvi = SMatrix{3, 3}(view(Mi, Nyi+1:Nyi+3, Nui+1:Nui+3))
        f_dωi = SMatrix{3, 3}(view(Mi, Nyi+1:Nyi+3, Nui+4:Nui+6))
        m_dλi = view(Mi, Nyi+4:Nyi+6, 1:Nui)
        m_dvi = SMatrix{3, 3}(view(Mi, Nyi+4:Nyi+6, Nui+1:Nui+3))
        m_dωi = SMatrix{3, 3}(view(Mi, Nyi+4:Nyi+6, Nui+4:Nui+6))
        # propgate derivatives using chain rule
        d_dV = d_dvi * dvi_dV
        d_dΩ = d_dvi * dvi_dΩ + d_dωi * dωi_dΩ
        f_dV = f_dvi * dvi_dV
        f_dΩ = f_dvi * dvi_dΩ + f_dωi * dωi_dΩ
        m_dV = m_dvi * dvi_dV
        m_dΩ = m_dvi * dvi_dΩ + m_dωi * dωi_dΩ
        # save aerodynamic input mass matrices
        M[iya[iyas[i]], iua[iuas[i]]] = d_dλi
        M[iya[iyas[i]], ius[7:9]] = d_dV
        M[iya[iyas[i]], ius[10:12]] = d_dΩ
        # section loads (in body frame)
        Fi_dλi = Cab*ΔL*f_dλi
        Fi_dV = Cab*ΔL*f_dV
        Fi_dΩ = Cab*ΔL*f_dΩ
        Mi_dλi = Cab*ΔL*m_dλi
        Mi_dV = Cab*ΔL*m_dV
        Mi_dΩ = Cab*ΔL*m_dΩ
        # add loads to total aerodynamic loads
        Fa_dλi = Fi_dλi
        Fa_dV += Fi_dV
        Fa_dΩ += Fi_dΩ
        Ma_dλi = GXBeam.tilde(pe)*Fi_dλi + Mi_dλi
        Ma_dV += GXBeam.tilde(pe)*Fi_dV + Mi_dV
        Ma_dΩ += GXBeam.tilde(pe)*Fi_dΩ + Mi_dΩ
        # save load mass matrix entries for this element
        M[iys[8:10], iua[iuas[i]]] = Fa_dλi
        M[iys[11:13], iua[iuas[i]]] = Ma_dλi
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

function get_inputs_from_state_rates(aero::LiftingLine{N,T}, stru::RigidBody,
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
    npc = number_of_parameters(aero, stru) # number of additional parameters for coupled model

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    ipc = npa + nps + 1 : npa + nps + npc # indices of additional parameters
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
    pc = view(p, ipc) # additional parameters for coupled model
    duas = view.(Ref(dua), iuas) # aerodynamic state rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states and rates
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = us
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dus

    # global parameters (first 18 additional parameters)
    Vx, Vy, Vz, ρ, g, m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz, Fx, Fy, Fz, Mx, My, Mz = pc

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - SVector(ur, vr, wr)
    dVinf = -SVector(dur, dvr, dwr)

    # body angular velocity
    Ω = SVector(pr, qr, rr)
    dΩ = SVector(dpr, dqr, drr)

    # aerodynamic forces and moments
    Fa = @SVector zeros(3)
    Ma = @SVector zeros(3)
    for i = 1:N
        # section aerodynamic model
        model = aero.models[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)
        # calculate section properties
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
        duai = SVector{Nui}(duas[i]) # aerodynamic rate variables
        uai = SVector{Nui}(uas[i]) # aerodynamic state variables
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters
        dvi = Cba*dVinf - Cba*cross(dΩ, pe) # local linear freestream acceleration
        dωi = Cba*dΩ # local angular freestream acceleration
        # calculate aerodynamic inputs and section loads
        dui = vcat(duai, dvi, dωi)
        ui = vcat(uai, vi, ωi) # section state variables
        pi = vcat(pai, ρ) # section parameters
        yi = get_inputs_from_state_rates(model, LiftingLineSection(),
            dui, ui, pi, t)
        # save aerodynamic inputs
        y[iya[iyas[i]]] .= view(yi, 1:Nyi)
        # section loads (in body frame)
        Fi = Cab*ΔL*SVector(yi[Nyi+1], yi[Nyi+2], yi[Nyi+3])
        Mi = Cab*ΔL*SVector(yi[Nyi+4], yi[Nyi+5], yi[Nyi+6])
        # add loads to total aerodynamic loads
        Fa += Fi
        Ma += cross(pe, Fi) + Mi
    end

    # save total forces and moments
    y[iys[8:10]] = Fa
    y[iys[11:13]] = Ma

    return y
end
