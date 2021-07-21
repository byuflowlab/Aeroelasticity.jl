"""
    couple_models(aero::LiftingLine, stru::GEBT, dyn::RigidBody)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a geometrically exact beam theory and rigid body model.  This model introduces
additional parameters corresponding to the freestream velocity components
``\\begin{bmatrix} V_x & V_y & V_z \\end{bmatrix}^T``, air density ``\\rho``,
gravitational constant ``g``, and external forces ``F_{x,i}``, ``F_{y,i}``,
``F_{z,i}``, ``M_{x,i}``, ``M_{y,i}``, ``M_{z,i}`` or displacements ``u_{x,i}``,
``u_{y,i}``, ``u_{z,i}``, ``\\theta_{x,i}``, ``\\theta_{y,i}``, ``\\theta_{z,i}``
applied to each node.

** When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface **
"""
couple_models(aero::LiftingLine, stru::GEBT, dyn::RigidBody) = (aero, stru, dyn)

# --- traits --- #

function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}, ::Type{<:RigidBody}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}, ::Type{<:RigidBody}) where {N,T}
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

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}, ::Type{<:RigidBody}) where {N,T}
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

function number_of_parameters(aero::LiftingLine, stru::GEBT, dyn::RigidBody)
    return 5 + 6*length(stru.icol_pt)
end

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, dyn::RigidBody, u, p, t) where {N,T}

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nud = number_of_states(dyn) # number of rigid body states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    npc = number_of_parameters(aero, stru, dyn) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iud = nua + nus + 1 : nua + nus + nud # indices of rigid body states
    iyd = nya + nys + 1 : nya + nys + nyd # indices of rigid body inputs
    ipd = npa + nps + 1 : npa + nps + npd # indices of rigid body parameters
    ipc = npa + nps + npd + 1 : npa + nps + npd + npc # indices of additional parameters
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
    ud = view(u, iud) # rigid body state variables
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    pc = view(p, ipc) # additional parameters
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = ud

    # global parameters (first 5 additional parameters)
    Vx, Vy, Vz, ρ, g = pc

    # body linear velocity
    V = SVector(ur, vr, wr)

    # body angular velocity
    Ω = SVector(pr, qr, rr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - V

    # gravity vector
    sϕ, cϕ = sincos(ϕr)
    sθ, cθ = sincos(θr)
    gvec = g*SVector(-sθ, cθ*sϕ, cθ*cϕ)

    # construct assembly from parameters
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # initialize rigid body properties
    mass = 0.0
    Ir = @SMatrix zeros(3,3)

    # save prescribed point loads/displacements
    for ip = 1:npoint
        poff = 5 + 6*(ip-1)
        yoff = 6*(ip-1)
        ys[yoff+1] = pc[poff+1]
        ys[yoff+2] = pc[poff+2]
        ys[yoff+3] = pc[poff+3]
        ys[yoff+4] = pc[poff+4]
        ys[yoff+5] = pc[poff+5]
        ys[yoff+6] = pc[poff+6]
    end

    # add prescribed loads to total forces and moments
    for ip = 1:npoint
        # point location
        point = assembly.points[ip]
        # prescribed condition identities
        isforce = stru.isforce[ip]
        # point loads
        poff = 5 + 6*(ip-1)
        F1 = ifelse(isforce[1], pc[poff+1], zero(pc[poff+1]))
        F2 = ifelse(isforce[2], pc[poff+2], zero(pc[poff+2]))
        F3 = ifelse(isforce[3], pc[poff+3], zero(pc[poff+3]))
        M1 = ifelse(isforce[4], pc[poff+4], zero(pc[poff+4]))
        M2 = ifelse(isforce[5], pc[poff+5], zero(pc[poff+5]))
        M3 = ifelse(isforce[6], pc[poff+6], zero(pc[poff+6]))
        Fp = SVector(F1, F2, F3)
        Mp = SVector(M1, M2, M3)
        Ftot += Fp
        Mtot += cross(point, Fp) + Mp
    end

    # calculate aerodynamic inputs, distributed loads, and element properties
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
        ΔL = element.L # element length
        pe = element.x + u_elem # location of element (in body frame)
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab # local to body transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1] # local structural to aerodynamic transformation
        vi = R*CtCab'*Vinf - R*GXBeam.element_linear_velocity(element, P_elem, H_elem) # local linear freestream velocity
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem) # local angular freestream velocity
        poff = 3*npoint + 36*(i-1)
        μ = ps[poff + 31] # mass per unit length
        xm2 = ps[poff + 32] # center of mass location
        xm3 = ps[poff + 33] # center of mass location
        i22 = ps[poff + 34] # rotational inertia
        i33 = ps[poff + 35] # rotational inertia
        i23 = ps[poff + 36] # product of inertia
        me = ΔL*μ # mass
        Ie = ΔL*(@SMatrix [i22+i33 0 0; 0 i22 -i23; 0 -i23 i33]) # local frame inertia
        Ie = (CtCab*R')*Ie*(R*CtCab') # body frame inertia
        Ie = me*(pe'*pe*I - pe*pe') + Ie # body frame inertia about origin
        λi = SVector{Nui}(uas[i]) # aerodynamic state variables
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters
        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        # calculate aerodynamic inputs and section loads
        ui = vcat(λi, vi, ωi) # section state variables
        pi = vcat(pai, ρ) # section parameters
        yi = get_inputs(model, LiftingLineSection(), ui, pi, t)

        # section aerodynamic inputs
        y[iya[iyas[i]]] .= view(yi, 1:Nyi)

        # section aerodynamic loads (in body frame)
        fi = CtCab*R'*SVector(yi[Nyi+1], yi[Nyi+2], yi[Nyi+3])
        mi = CtCab*R'*SVector(yi[Nyi+4], yi[Nyi+5], yi[Nyi+6])

        # add gravitational distributed loads (in body frame)
        fi += μ*gvec
        mi += cross(CtCab*R'*SVector(0, xm2, xm3), μ*gvec)

        # save distributed loads for this element (in body frame)
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1 : yoff+3]] = fi
        y[iys[yoff+4 : yoff+6]] = mi

        # add distributed loads to total forces and moments
        Ftot += ΔL*fi
        Mtot += cross(pe, ΔL*fi) + ΔL*mi

        # add element mass and inertia to total mass and inertia
        mass += me
        Ir += Ie
    end

    # save body frame linear and angular velocities
    yoff = 6*npoint + 6*nelem
    y[iys[yoff+1:yoff+3]] = V
    y[iys[yoff+4:yoff+6]] = Ω

    # save rigid body model inputs
    y[iyd[1]] = mass
    y[iyd[2]] = Ixx = Ir[1,1]
    y[iyd[3]] = Iyy = Ir[2,2]
    y[iyd[4]] = Izz = Ir[3,3]
    y[iyd[5]] = Ixz = Ir[1,3]
    y[iyd[6]] = Ixy = Ir[1,2]
    y[iyd[7]] = Iyz = Ir[2,3]
    y[iyd[8:10]] = Ftot
    y[iyd[11:13]] = Mtot

    return y
end

function get_input_mass_matrix!(My, aero::LiftingLine{N,T}, stru::GEBT,
    dyn::RigidBody, u, p, t) where {N,T}

    # start with zero valued mass matrix
    My .= 0

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nud = number_of_states(dyn) # number of rigid body states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    npc = number_of_parameters(aero, stru, dyn) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iud = nua + nus + 1 : nua + nus + nud # indices of rigid body states
    iyd = nya + nys + 1 : nya + nys + nyd # indices of rigid body inputs
    ipd = npa + nps + 1 : npa + nps + npd # indices of rigid body parameters
    ipc = npa + nps + npd + 1 : npa + nps + npd + npc # indices of additional parameters
    iuas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    ua = view(u, iua) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    us = view(u, ius) # structural state variables
    ps = view(p, ips) # structural parameters
    ud = view(u, iud) # rigid body state variables
    pd = view(p, ipd) # rigid body parameters
    pc = view(p, ipc) # additional parameters
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = ud

    # global parameters (first 5 additional parameters)
    Vx, Vy, Vz, ρ, g = pc

    # body angular velocity
    Ω = SVector(pr, qr, rr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - SVector(ur, vr, wr)

    # freestream acceleration
    dVinf_dV = -I

    # construct assembly from parameters
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # initialize total forces and moments
    Ftot_dV = @SMatrix zeros(3, 3)
    Ftot_dΩ = @SMatrix zeros(3, 3)
    Mtot_dV = @SMatrix zeros(3, 3)
    Mtot_dΩ = @SMatrix zeros(3, 3)

    # calculate aerodynamic inputs, distributed loads, and element properties
    for i = 1:N
        # aerodynamic model and structural element
        model = aero.models[i]
        element = assembly.elements[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)

        # structural state variables
        icol = stru.icol_beam[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = GXBeam.rotation_parameter_scaling(θ_elem)
        θ_elem *= scaling # angular displacement (Wiener-Milenkovic parameters)

        # length and location of the element
        ΔL = element.L
        pe = element.x + u_elem # location of element (in body frame)

        # local to body transformation
        CtCab = GXBeam.get_C(θ_elem)'*element.Cab

        # local structural to aerodynamic transformation
        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # local freestream linear and angular velocities
        vi = R*CtCab'*Vinf - R*GXBeam.element_linear_velocity(element, P_elem, H_elem)
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem)

        # local freestream linear and angular accelerations
        dvi_dV = -R*CtCab'
        dvi_dΩ = R*CtCab'*GXBeam.tilde(pe)
        dvi_dPi = -R * element.minv11 * stru.mass_scaling
        dvi_dHi = -R * element.minv12 * stru.mass_scaling

        dωi_dΩ = R*CtCab'
        dωi_dPi = R * element.minv12' * stru.mass_scaling
        dωi_dHi = R * element.minv22 * stru.mass_scaling

        # aerodynamic states and parameters
        λi = SVector{Nui}(uas[i]) # aerodynamic state variables
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters

        # element inertial properties
        poff = 3*npoint + 36*(i-1)
        μ = ps[poff + 31] # mass per unit length
        xm2 = ps[poff + 32] # center of mass location
        xm3 = ps[poff + 33] # center of mass location
        i22 = ps[poff + 34] # rotational inertia
        i33 = ps[poff + 35] # rotational inertia
        i23 = ps[poff + 36] # product of inertia
        me = ΔL*μ # mass
        Ie = ΔL*(@SMatrix [i22+i33 0 0; 0 i22 -i23; 0 -i23 i33]) # local frame inertia
        Ie = (CtCab*R')*Ie*(R*CtCab') # body frame inertia
        Ie = me*(pe'*pe*I - pe*pe') + Ie # body frame inertia about origin

        # calculate lifting line section input mass matrix
        ui = vcat(λi, vi, ωi) # section state variables
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
        d_dV = d_dvi * dvi_dV
        d_dΩ = d_dvi * dvi_dΩ + d_dωi * dωi_dΩ
        d_dPi = d_dvi * dvi_dPi + d_dωi * dωi_dPi
        d_dHi = d_dvi * dvi_dHi + d_dωi * dωi_dHi
        f_dV = f_dvi * dvi_dV
        f_dΩ = f_dvi * dvi_dΩ + f_dωi * dωi_dΩ
        f_dPi = f_dvi * dvi_dPi + f_dωi * dωi_dPi
        f_dHi = f_dvi * dvi_dHi + f_dωi * dωi_dHi
        m_dV = m_dvi * dvi_dV
        m_dΩ = m_dvi * dvi_dΩ + m_dωi * dωi_dΩ
        m_dPi = m_dvi * dvi_dPi + m_dωi * dωi_dPi
        m_dHi = m_dvi * dvi_dHi + m_dωi * dωi_dHi

        # save aerodynamic input mass matrix entries
        icol = stru.icol_beam[i]
        My[iya[iyas[i]], iua[iuas[i]]] = d_dλi
        My[iya[iyas[i]], ius[icol+12:icol+14]] = d_dPi
        My[iya[iyas[i]], ius[icol+15:icol+17]] = d_dHi
        My[iya[iyas[i]], iud[7:9]] = d_dV
        My[iya[iyas[i]], iud[10:12]] = d_dΩ

        # section aerodynamic load mass matrices
        f_dλi = CtCab*R'*f_dλi
        f_dPi = CtCab*R'*f_dPi
        f_dHi = CtCab*R'*f_dHi
        f_dV = CtCab*R'*f_dV
        f_dΩ = CtCab*R'*f_dΩ
        m_dλi = CtCab*R'*m_dλi
        m_dPi = CtCab*R'*m_dPi
        m_dHi = CtCab*R'*m_dHi
        m_dV = CtCab*R'*m_dV
        m_dΩ = CtCab*R'*m_dΩ

        # # add apparent forces due to body frame linear and angular acceleration
        f_dV -= me*I
        f_dΩ -= me*GXBeam.tilde(pe)
        m_dΩ -= me*I

        # save load mass matrix entries
        offset = 6*npoint + 6*(i-1)
        My[iys[offset+1 : offset+3], iua[iuas[i]]] = f_dλi
        My[iys[offset+4 : offset+6], iua[iuas[i]]] = m_dλi
        My[iys[offset+1 : offset+3], ius[icol+12:icol+14]] = f_dPi
        My[iys[offset+1 : offset+3], ius[icol+15:icol+17]] = f_dHi
        My[iys[offset+4 : offset+6], ius[icol+12:icol+14]] = m_dPi
        My[iys[offset+4 : offset+6], ius[icol+15:icol+17]] = m_dHi
        My[iys[offset+1 : offset+3], iud[7:9]] = f_dV
        My[iys[offset+1 : offset+3], iud[10:12]] = f_dΩ
        My[iys[offset+4 : offset+6], iud[7:9]] = m_dV
        My[iys[offset+4 : offset+6], iud[10:12]] = m_dΩ

        # add contributions to total forces and moments
        My[iyd[8:10], iua[iuas[i]]] = ΔL*f_dλi
        My[iyd[8:10], ius[icol+12:icol+14]] = ΔL*f_dPi
        My[iyd[8:10], ius[icol+15:icol+17]] = ΔL*f_dHi
        Ftot_dV += ΔL*f_dV
        Ftot_dΩ += ΔL*f_dΩ
        My[iyd[11:13], iua[iuas[i]]] = GXBeam.tilde(pe)*ΔL*f_dλi + ΔL*m_dλi
        My[iyd[11:13], ius[icol+12:icol+14]] = GXBeam.tilde(pe)*ΔL*f_dPi + ΔL*m_dPi
        My[iyd[11:13], ius[icol+15:icol+17]] = GXBeam.tilde(pe)*ΔL*f_dHi + ΔL*m_dHi
        Mtot_dV += GXBeam.tilde(pe)*ΔL*f_dV + ΔL*m_dV
        Mtot_dΩ += GXBeam.tilde(pe)*ΔL*f_dΩ + ΔL*m_dΩ
    end

    # save rigid body model inputs
    My[iyd[8:10], iud[7:9]] = Ftot_dV
    My[iyd[8:10], iud[10:12]] = Ftot_dΩ
    My[iyd[11:13], iud[7:9]] = Mtot_dV
    My[iyd[11:13], iud[10:12]] = Mtot_dΩ

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{N,T}, stru::GEBT,
    dyn::RigidBody, du, u, p, t) where {N,T}

    # initialize input vector
    models = (aero, stru, dyn)
    TF = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Ny = number_of_inputs(models)
    y = zeros(TF, Ny)

    # extract number of state variables, inputs, and parameters
    nua = number_of_states(aero) # number of aerodynamic states
    nus = number_of_states(stru) # number of structural states
    nud = number_of_states(dyn) # number of rigid body states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    npc = number_of_parameters(aero, stru, dyn) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    iua = 1:nua # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ius = nua + 1 : nua + nus # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iud = nua + nus + 1 : nua + nus + nud # indices of rigid body states
    iyd = nya + nys + 1 : nya + nys + nyd # indices of rigid body inputs
    ipd = npa + nps + 1 : npa + nps + npd # indices of rigid body parameters
    ipc = npa + nps + npd + 1 : npa + nps + npd + npc # indices of additional parameters
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
    dud = view(du, iud) # rigid body state rates
    ud = view(u, iud) # rigid body state variables
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    pc = view(p, ipc) # additional parameters
    duas = view.(Ref(dua), iuas) # aerodynamic state rates for each section
    uas = view.(Ref(ua), iuas) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = ud

    # rigid body state rates
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dud

    # global parameters (first 5 additional parameters)
    Vx, Vy, Vz, ρ, g = pc

    # body linear velocity
    V = SVector(ur, vr, wr)

    # body linear acceleration
    dV = SVector(dur, dvr, dwr)

    # body angular velocity
    Ω = SVector(pr, qr, rr)

    # body angular acceleration
    dΩ = SVector(dpr, dqr, drr)

    # freestream velocity
    Vinf = SVector(Vx, Vy, Vz) - V

    # freestream acceleration
    dVinf = -dV

    # construct assembly from parameters
    npoint = length(stru.icol_pt)
    nelem = length(stru.icol_beam)
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # calculate aerodynamic inputs and distributed loads
    for i = 1:N
        # aerodynamic model and structural element
        model = aero.models[i]
        element = assembly.elements[i]
        Nui = number_of_states(model)
        Nyi = number_of_inputs(model)
        Npi = number_of_parameters(model)

        # structural state variables
        icol = stru.icol_beam[i]
        u_elem = SVector(us[icol], us[icol+1], us[icol+2]) # linear displacement
        θ_elem = SVector(us[icol+3], us[icol+4], us[icol+5]) # angular displacement
        F_elem = SVector(us[icol+6], us[icol+7], us[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(us[icol+9], us[icol+10], us[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(us[icol+12], us[icol+13], us[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(us[icol+15], us[icol+16], us[icol+17]) .* stru.mass_scaling # angular momentum

        dF_elem = SVector(dus[icol+6], dus[icol+7], dus[icol+8]) .* stru.force_scaling
        dM_elem = SVector(dus[icol+9], dus[icol+10], dus[icol+11]) .* stru.force_scaling
        dP_elem = SVector(dus[icol+12], dus[icol+13], dus[icol+14]) .* stru.mass_scaling
        dH_elem = SVector(dus[icol+15], dus[icol+16], dus[icol+17]) .* stru.mass_scaling

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
        vi = R*CtCab'*Vinf - R*GXBeam.element_linear_velocity(element, P_elem, H_elem) # local linear freestream velocity
        ωi = R*GXBeam.element_angular_velocity(element, P_elem, H_elem) # local angular freestream velocity

        # local freestream linear and angular accelerations
        dvi = R*CtCab'*dVinf - R*CtCab'*cross(dΩ, pe) - R*GXBeam.element_linear_velocity(element, dP_elem, dH_elem)
        dωi = R*CtCab'*dΩ + R*GXBeam.element_angular_velocity(element, dP_elem, dH_elem)

        # aerodynamic rates, states, and parameters
        dλi = SVector{Nui}(duas[i]) # aerodynamic rates
        λi = SVector{Nui}(uas[i]) # aerodynamic states
        pai = SVector{Npi}(pas[i]) # aerodynamic parameters

        # element inertial properties
        poff = 3*npoint + 36*(i-1)
        μ = ps[poff + 31] # mass per unit length
        xm2 = ps[poff + 32] # center of mass location
        xm3 = ps[poff + 33] # center of mass location
        i22 = ps[poff + 34] # rotational inertia
        i33 = ps[poff + 35] # rotational inertia
        i23 = ps[poff + 36] # product of inertia
        me = ΔL*μ # element mass
        Ie = ΔL*(@SMatrix [i22+i33 0 0; 0 i22 -i23; 0 -i23 i33]) # local frame inertia
        Ie = (CtCab*R')*Ie*(R*CtCab') # body frame inertia
        Ie = me*(pe'*pe*I - pe*pe') + Ie # body frame inertia about origin

        # calculate lifting line section inputs from state rates
        dui = vcat(dλi, dvi, dωi) # section state rates
        ui = vcat(λi, vi, ωi) # section states
        pi = vcat(pai, ρ) # section parameters
        yi = get_inputs_from_state_rates(model, LiftingLineSection(),
            dui, ui, pi, t)

        # section aerodynamic inputs
        y[iya[iyas[i]]] .= view(yi, 1:Nyi)

        # section aerodynamic loads
        fi = CtCab*R'*SVector(yi[Nyi+1], yi[Nyi+2], yi[Nyi+3])
        mi = CtCab*R'*SVector(yi[Nyi+4], yi[Nyi+5], yi[Nyi+6])

        # add apparent forces due to body frame linear and angular acceleration
        fi += me*dV - cross(me*dΩ, pe)
        mi += me*dΩ

        # section distributed loads
        yoff = 6*npoint + 6*(i-1)
        y[iys[yoff+1:yoff+3]] = fi
        y[iys[yoff+4:yoff+6]] = mi

        # add distributed loads to total forces and moments
        Ftot += ΔL*fi
        Mtot += cross(pe, ΔL*fi) + ΔL*mi
    end

    # save state rate contribution to body frame linear and angular velocities
    ys[6*N+1 : 6*N+3] .= 0
    ys[6*N+4 : 6*N+6] .= 0

    # save state rate contribution to rigid body model inputs
    yd[1] = 0
    yd[2] = Ixx = 0
    yd[3] = Iyy = 0
    yd[4] = Izz = 0
    yd[5] = Ixz = 0
    yd[6] = Ixy = 0
    yd[7] = Iyz = 0
    yd[8:10] = Ftot
    yd[11:13] = Mtot

    return y
end
