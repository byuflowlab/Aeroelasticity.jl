"""
    couple_models(aero::LiftingLine, stru::GEBT)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a geometrically exact beam theory model.  This model introduces additional
parameters corresponding to the freestream velocity components ``V_x``, ``V_y``,
``V_z``, air density ``\\rho``, body linear velocity components ``u``,
``v``, ``w``, body angular velocity components ``p``, ``q``, ``r``, external
loads ``F_{x,i}``, ``F_{y,i}``, ``F_{z,i}``, ``M_{x,i}``, ``M_{y,i}``,
``M_{z,i}`` or displacements ``u_{x,i}``, ``u_{y,i}``, ``u_{z,i}``,
``\\theta_{x,i}``, ``\\theta_{y,i}``, ``\\theta_{z,i}`` for each node, and
constant distributed loads ``f_{x,i}``, ``f_{y,i}``, ``f_{z,i}``, ``m_{x,i}``,
``m_{y,i}``, ``m_{z,i}`` for each beam element (excluding aerodynamic loads).

** When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface **
"""
couple_models(aero::LiftingLine, stru::GEBT)

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
    return 10 + 6*length(stru.icol_pt) + 6*length(stru.icol_beam)
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
    npc = number_of_parameters(aero, stru) # number of additional parameters

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

    # global parameters (first 10 additional parameters)
    Vx, Vy, Vz, ρ, ur, vr, wr, pr, qr, rr = pc

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - SVector(ur, vr, wr)

    # body angular velocity
    Ω = SVector(pr, qr, rr)

    # construct assembly from structural parameters
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # save prescribed point loads/displacements
    for ip = 1:npoint
        yoff = nya + 6*(ip-1)
        poff = npa + nps + 10 + 6*(ip-1)
        y[yoff+1] = p[poff+1]
        y[yoff+2] = p[poff+2]
        y[yoff+3] = p[poff+3]
        y[yoff+4] = p[poff+4]
        y[yoff+5] = p[poff+5]
        y[yoff+6] = p[poff+6]
    end

    # save aerodynamic inputs and distributed loads
    for i = 1:N
        # aerodynamic model and structural element
        model = aero.models[i]
        element = assembly.elements[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)
        # section properties
        icol = stru.icol_beam[i]
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
        uai = SVector{Nui}(uas[i]) # aerodynamic state variables
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters
        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.
        # calculate aerodynamic inputs and section loads
        ui = vcat(uai, vi, ωi) # section state variables
        pi = vcat(pai, ρ) # section parameters
        yi = get_inputs(model, LiftingLineSection(), ui, pi, t)
        # section aerodynamic inputs
        y[iya[iyas[i]]] .= view(yi, 1:Nyi)
        # section aerodynamic loads (in body frame)
        fi = CtCab*R'*SVector(yi[Nyi+1], yi[Nyi+2], yi[Nyi+3])
        mi = CtCab*R'*SVector(yi[Nyi+4], yi[Nyi+5], yi[Nyi+6])
        # add constant distributed loads (in body frame)
        poff = 10 + 6*npoint + 6*(i-1)
        fi += SVector(pc[poff+1], pc[poff+2], pc[poff+3])
        mi += SVector(pc[poff+4], pc[poff+5], pc[poff+6])
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
    pa = view(p, ipa) # aerodynamic parameters
    us = view(u, ius) # structural state variables
    ps = view(p, ips) # structural parameters
    pc = view(p, ipc) # additional parameters
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # global parameters (first 10 additional parameters)
    Vx, Vy, Vz, ρ, ur, vr, wr, pr, qr, rr = pc

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
        # aerodynamic model and structural element
        model = aero.models[i]
        element = assembly.elements[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)
        # section properties
        icol = stru.icol_beam[i]
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
        uai = SVector{Nui}(uas[i]) # aerodynamic state variables
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters
        dvi_dPi = -R * element.minv11 * stru.mass_scaling
        dvi_dHi = -R * element.minv12 * stru.mass_scaling
        dωi_dPi = R * element.minv12' * stru.mass_scaling
        dωi_dHi = R * element.minv22 * stru.mass_scaling
        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.
        # calculate aerodynamic input and section loads mass matrices
        ui = vcat(uai, vi, ωi) # section state variables
        pi = vcat(pai, ρ) # section parameters
        Myi = get_input_mass_matrix(model, LiftingLineSection(), ui, pi, t)
        # separate into component mass matrices
        d_dλi = SMatrix{Nyi,Nui}(view(Myi, 1:Nyi, 1:Nui))
        d_dvi = SMatrix{Nyi,3}(view(Myi, 1:Nyi, Nui+1:Nui+3))
        d_dωi = SMatrix{Nyi,3}(view(Myi, 1:Nyi, Nui+4:Nui+6))
        f_dλi = view(Myi, Nyi+1:Nyi+3, 1:Nui)
        f_dvi = SMatrix{3, 3}(view(Myi, Nyi+1:Nyi+3, Nui+1:Nui+3))
        f_dωi = SMatrix{3, 3}(view(Myi, Nyi+1:Nyi+3, Nui+4:Nui+6))
        m_dλi = view(Myi, Nyi+4:Nyi+6, 1:Nui)
        m_dvi = SMatrix{3, 3}(view(Myi, Nyi+4:Nyi+6, Nui+1:Nui+3))
        m_dωi = SMatrix{3, 3}(view(Myi, Nyi+4:Nyi+6, Nui+4:Nui+6))
        # propagate derivatives using chain rule
        d_dPi = d_dvi * dvi_dPi + d_dωi * dωi_dPi
        d_dHi = d_dvi * dvi_dHi + d_dωi * dωi_dHi
        f_dPi = f_dvi * dvi_dPi + f_dωi * dωi_dPi
        f_dHi = f_dvi * dvi_dHi + f_dωi * dωi_dHi
        m_dPi = m_dvi * dvi_dPi + m_dωi * dωi_dPi
        m_dHi = m_dvi * dvi_dHi + m_dωi * dωi_dHi
        # save aerodynamic input mass matrix entries
        icol = stru.icol_beam[i]
        My[iya[iyas[i]], iua[iuas[i]]] = d_dλi
        My[iya[iyas[i]], ius[icol+12]:ius[icol+14]] = d_dPi
        My[iya[iyas[i]], ius[icol+15]:ius[icol+17]] = d_dHi
        # save load mass matrix entries (in body frame)
        offset = 6*npoint + 6*(i-1)
        My[iys[offset+1 : offset+3], iua[iuas[i]]] = CtCab*R'*f_dλi
        My[iys[offset+4 : offset+6], iua[iuas[i]]] = CtCab*R'*m_dλi
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

    # global parameters (first 10 additional parameters)
    Vx, Vy, Vz, ρ, ur, vr, wr, pr, qr, rr = pc

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
        # aerodynamic model and structural element
        model = aero.models[i]
        element = assembly.elements[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)

        # section properties
        icol = stru.icol_beam[i]
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

        # aerodynamic rates, states, and parameters
        duai = SVector{Nui}(duas[i]) # aerodynamic rates
        uai = SVector{Nui}(uas[i]) # aerodynamic states
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters


        # calculate lifting line section inputs from state rates
        dui = vcat(duai, dvi, dωi) # section state rates
        ui = vcat(uai, vi, ωi) # section states
        pi = vcat(pai, ρ) # section parameters
        yi = get_inputs_from_state_rates(model, LiftingLineSection(),
            dui, ui, pi, t)

        # section aerodynamic inputs
        y[iya[iyas[i]]] .= view(yi, 1:Nyi)

        # section aerodynamic loads
        fi = CtCab*R'*SVector(yi[Nyi+1], yi[Nyi+2], yi[Nyi+3])
        mi = CtCab*R'*SVector(yi[Nyi+4], yi[Nyi+5], yi[Nyi+6])

        # section distributed loads
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1:yoff+3]] = fi
        y[iys[yoff+4:yoff+6]] = mi
    end

    return y
end
