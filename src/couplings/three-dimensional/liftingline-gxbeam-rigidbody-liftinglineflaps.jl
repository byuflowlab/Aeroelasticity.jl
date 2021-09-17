"""
    couple_models(aero::LiftingLine, stru::GEBT, dyn::RigidBody, flap::LiftingLineFlaps)

Create a coupled model using a lifting line aerodynamic model, a
geometrically exact beam theory model, a rigid body dynamics model, and a lifting
line control surface model.  This model introduces additional parameters
corresponding to the freestream air density ``\\rho_\\infty``, followed by the
external loads ``F_{x,i}, F_{y,i}, F_{z,i}, M_{x,i}, M_{y,i}, M_{z,i}`` or
displacements ``u_{x,i}, u_{y,i}, u_{z,i}, \\theta_{x,i}, \\theta_{y,i},
\\theta_{z,i}`` for each node, followed by the constant distributed loads
``f_{x,i}, f_{y,i}, f_{z,i}, m_{x,i}, m_{y,i}, m_{z,i}`` applied on each beam
element (excluding aerodynamic loads), followed by the control surface
deflections ``\\delta_1, \\delta_2, \\dots, \\delta_N``.

**NOTE: When using this model, the local frame for each beam element should be
oriented with the x-axis along the beam's axis, the y-axis forward, and the
z-axis normal to the surface**
"""
function couple_models(aero::LiftingLine, stru::GEBT, dyn::RigidBody,
    flap::LiftingLineFlaps)

    return (aero, stru, dyn, flap)
end

# --- traits --- #

function number_of_additional_parameters(aero::LiftingLine, stru::GEBT, dyn::RigidBody,
    flap::LiftingLineFlaps)

    return 1 + 6*length(stru.icol_point) + 6*length(stru.icol_elem) + length(flap.gains)
end

function coupling_inplaceness(::Type{<:LiftingLine}, ::Type{<:GEBT}, ::Type{<:RigidBody},
    ::Type{<:LiftingLineFlaps})

    return InPlace()
end

function coupling_rate_jacobian_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT},
    ::Type{<:RigidBody}, ::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NA,NF,NG,TA,TF,TG}

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
    ::Type{<:RigidBody}, ::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NA,NF,NG,TA,TF,TG}

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
    ::Type{RigidBody}, ::Type{LiftingLineFlaps}) = Nonlinear()


function coupling_time_gradient_type(::Type{LiftingLine{NA,TA}}, ::Type{<:GEBT},
    ::Type{<:RigidBody}, ::Type{LiftingLineFlaps{NF,NG,TF,TG}}) where {NA,NF,NG,TA,TF,TG}

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

function get_coupling_inputs!(y, aero::LiftingLine{NA,TA}, stru::GEBT, dyn::RigidBody,
    flap::LiftingLineFlaps{NF,NG,TF,TG}, dx, x, p, t) where {NA,NF,NG,TA,TF,TG}

    # extract number of state variables, inputs, and parameters
    nxa = number_of_states(aero) # number of aerodynamic states
    nxs = number_of_states(stru) # number of structural states
    nxd = number_of_states(dyn) # number of rigid body states
    nxf = number_of_states(flap) # number of control surface states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    nyd = number_of_inputs(dyn) # number of rigid body inputs
    nyf = number_of_inputs(flap) # number of control surface inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    npd = number_of_parameters(dyn) # number of rigid body parameters
    npf = number_of_parameters(flap) # number of control surface parameters
    npadd = number_of_additional_parameters(aero, stru, dyn, flap) # number of additional parameters

    # get indices for state variables, inputs, and parameters
    ixa = 1:nxa # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    ixs = nxa + 1 : nxa + nxs # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    ixd = nxa + nxs + 1 : nxa + nxs + nxd # indices of rigid body states
    iyd = nya + nys + 1 : nya + nys + nyd # indices of rigid body inputs
    ipd = npa + nps + 1 : npa + nps + npd # indices of rigid body parameters
    ixf = nxa + nxs + nxd + 1 : nxa + nxs + nxd + nxf # indices of control surface states
    iyf = nya + nys + nyd + 1 : nya + nys + nyd + nyf # indices of control surface inputs
    ipf = npa + nps + npd + 1 : npa + nps + npd + npf # indices of control surface parameters
    ipadd = npa + nps + npd + npf + 1 : npa + nps + npd + npf + npadd # indices of additional parameters
    ixas = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section
    ixfs = state_indices(flap.models) # indices of control surface states for each section
    iyfs = input_indices(flap.models) # indices of control surface inputs for each section
    ipfs = parameter_indices(flap.models) # indices of control surface parameters for each section

    # separate state variables, inputs, and parameters
    dxa = view(dx, ixa) # aerodynamic rates
    xa = view(x, ixa) # aerodynamic states
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dxs = view(dx, ixs) # structural rates
    xs = view(x, ixs) # structural states
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    dxd = view(dx, ixd) # rigid body rates
    xd = view(x, ixd) # rigid body states
    yd = view(y, iyd) # rigid body inputs
    pd = view(p, ipd) # rigid body parameters
    dxf = view(dx, ixf) # structural rates
    xf = view(x, ixf) # structural states
    yf = view(y, iyf) # structural inputs
    pf = view(p, ipf) # structural parameters
    padd = view(p, ipadd) # additional parameters for coupled model
    dxas = view.(Ref(dxa), ixas) # aerodynamic rates for each section
    xas = view.(Ref(xa), ixas) # aerodynamic states for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section
    dxfs = view.(Ref(dxf), ixfs) # control surface rates for each section
    xfs = view.(Ref(xf), ixfs) # control surface states for each section
    yfs = view.(Ref(yf), iyfs) # control surface inputs for each section
    pfs = view.(Ref(pf), ipfs) # control surface parameters for each section

    # number of points and elements
    npoint = length(stru.icol_point)
    nelem = length(stru.icol_elem)

    # rigid body states
    xr, yr, zr, ϕr, θr, ψr, ur, vr, wr, pr, qr, rr = xd

    # rigid body state rates
    dxr, dyr, dzr, dϕr, dθr, dψr, dur, dvr, dwr, dpr, dqr, drr = dxd

    # extract global parameters
    ρ = padd[1] # freestream air density
    δ = padd[SVector{NG}(1 + 6*npoint + 6*nelem + 1 : 1 + 6*npoint + 6*nelem + NG)]
    dδ = zero(δ)

    # rigid body linear and angular velocity
    V = SVector(ur, vr, wr)
    Ω = SVector(pr, qr, rr)

    # rigid body linear and angular acceleration
    dV = SVector(dur, dvr, dwr)
    dΩ = SVector(dpr, dqr, drr)

    # construct assembly from parameters
    assembly = gxbeam_assembly(ps, npoint, nelem, stru.start, stru.stop)

    # initialize total forces and moments
    Ftot = @SVector zeros(3)
    Mtot = @SVector zeros(3)

    # initialize rigid body properties
    mass = 0.0
    Ir = @SMatrix zeros(3,3)

    # save prescribed point loads/displacements
    for ip = 1:npoint
        yoff = 6*(ip-1)
        poff = 1 + 6*(ip-1)
        ys[yoff+1] = padd[poff+1]
        ys[yoff+2] = padd[poff+2]
        ys[yoff+3] = padd[poff+3]
        ys[yoff+4] = padd[poff+4]
        ys[yoff+5] = padd[poff+5]
        ys[yoff+6] = padd[poff+6]
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
        Nxai = number_of_states(section_aero)
        Nyai = number_of_inputs(section_aero)
        Npai = number_of_parameters(section_aero)
        Nxsi = number_of_states(section_stru)
        Nysi = number_of_inputs(section_stru)
        Npsi = number_of_parameters(section_stru)
        Nxfi = number_of_states(section_flap)
        Nyfi = number_of_inputs(section_flap)
        Npfi = number_of_parameters(section_flap)
        Nxci = number_of_states(section_ctrl)
        Nyci = number_of_inputs(section_ctrl)
        Npci = number_of_parameters(section_ctrl)

        # section properties
        icol = stru.icol_elem[i]
        element = assembly.elements[i]
        u_elem = SVector(xs[icol], xs[icol+1], xs[icol+2]) # linear displacement
        θ_elem = SVector(xs[icol+3], xs[icol+4], xs[icol+5]) # angular displacement
        F_elem = SVector(xs[icol+6], xs[icol+7], xs[icol+8]) .* stru.force_scaling # internal force
        M_elem = SVector(xs[icol+9], xs[icol+10], xs[icol+11]) .* stru.force_scaling # internal moment
        P_elem = SVector(xs[icol+12], xs[icol+13], xs[icol+14]) .* stru.mass_scaling # linear momentum
        H_elem = SVector(xs[icol+15], xs[icol+16], xs[icol+17]) .* stru.mass_scaling # angular momentum

        dF_elem = SVector(dxs[icol+6], dxs[icol+7], dxs[icol+8]) .* stru.force_scaling
        dM_elem = SVector(dxs[icol+9], dxs[icol+10], dxs[icol+11]) .* stru.force_scaling
        dP_elem = SVector(dxs[icol+12], dxs[icol+13], dxs[icol+14]) .* stru.mass_scaling
        dH_elem = SVector(dxs[icol+15], dxs[icol+16], dxs[icol+17]) .* stru.mass_scaling

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
        dvi = -R*(CtCab'*dV + CtCab'*cross(dΩ, pe) + GXBeam.element_linear_velocity(element, dP_elem, dH_elem))
        dωi = R*CtCab'*dΩ + R*GXBeam.element_angular_velocity(element, dP_elem, dH_elem)

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

        # NOTE: We assume that the local beam frame y-axis is oriented towards
        # the leading edge and the z-axis is oriented up
        # NOTE: We assume the local aerodynamic frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        # section state variable rates
        duai = SVector{Nxai}(dxas[i]) # aerodynamic state variables
        dusi = vcat(dvi, dωi) # structural state variables
        dufi = SVector{Nxfi}(dxfs[i]) # control surface state variables
        duci = SMatrix{1,NG}(getindex.(flap.gains, i)) * dδ # controller state variables
        dui = vcat(duai, dusi, dufi, duci)

        # section state variables
        uai = SVector{Nxai}(xas[i]) # aerodynamic state variables
        usi = vcat(vi, ωi) # structural state variables
        ufi = SVector{Nxfi}(xfs[i]) # control surface state variables
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
        fi = CtCab*R'*SVector(ysi[1], ysi[2], ysi[3])
        mi = CtCab*R'*SVector(ysi[4], ysi[5], ysi[6])

        # add apparent forces due to body frame linear and angular acceleration
        fi += me*dV - cross(me*dΩ, pe)
        mi += me*dΩ

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

        # add element mass and inertia to total mass and inertia
        mass += me
        Ir += Ie
    end

    # save body frame linear and angular velocities
    yoff = 6*npoint + 6*nelem
    y[iys[yoff+1]] = ur
    y[iys[yoff+2]] = vr
    y[iys[yoff+3]] = wr
    y[iys[yoff+4]] = pr
    y[iys[yoff+5]] = qr
    y[iys[yoff+6]] = rr

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

# --- performance overloads --- #

# TODO


# --- convenience methods --- #

function set_additional_parameters!(padd, model::LiftingLine, stru::GEBT,
    dyn::RigidBody, flap::LiftingLineFlaps{NF,NG,TF,TG}; rho, point_conditions,
    element_loads, delta) where {NF,NG,TF,TG}

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

    padd[1 + 6*np + 6*ne + 1 : 1 + 6*np + 6*ne + NG] = delta

    return padd
end

function separate_additional_parameters(model::LiftingLine, stru::GEBT,
    dyn::RigidBody, flap::LiftingLineFlaps{NF,NG,TF,TG}, padd) where {NF,NG,TF,TG}

    np = length(stru.icol_point)
    ne = length(stru.icol_elem)

    rho = padd[1]

    point_conditions = reshape(view(padd, 1 + 1 : 1 + 6*np), 6, np)

    element_loads = reshape(view(padd, 1 + 6*np + 1 : 1 + 6*np + 6*ne), 6, ne)

    delta = view(padd, 1 + 6*np + 6*ne + 1 : 1 + 6*np + 6*ne + NG)

    return (rho = rho, point_conditions = point_conditions, element_loads = element_loads,
        delta = delta)
end
