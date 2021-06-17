"""
    LiftingLine{N,T} <: AbstractModel

Lifting line model with `N` cross sections, using the aerodynamic models in `T`.
State variables, inputs, and parameters correspond to the state variables, inputs,
and parameters of each of the cross sections concatenated
"""
struct LiftingLine{N,T} <: AbstractModel
    models::T
end

"""
    LiftingLineSection <: AbstractModel

Lifting line section model with state variables ``q = \\begin{bmatrix} v_x & v_y &
v_z & \\omega_x & \\omega_y & \\omega_z \\end{bmatrix}^T``, inputs
``r = \\begin{bmatrix} L' & Y' & D' & Mx, My, Mz \\end{bmatrix}^T``, and no
parameters.  Two-dimensional aerodynamic models may be extended to
three-dimensional models by coupling with this model.  Note that
this model has no rate equations of its own. Its state variables are
defined as functions of a 3D structural model's state variables.
"""
struct LiftingLineSection <: AbstractModel end

number_of_states(::Type{<:LiftingLineSection}) = 6
number_of_inputs(::Type{<:LiftingLineSection}) = 6
number_of_parameters(::Type{<:LiftingLineSection}) = 0
inplaceness(::Type{LiftingLineSection}) = OutOfPlace()

# --- Constructors --- #

"""
    LiftingLine(models)

Construct a lifting line aerodynamic model given a tuple of aerodynamic models.
"""
LiftingLine(models::T) where T<:NTuple{N,Any} where N = LiftingLine{N,T}(models)

LiftingLine{N}(models::T) where T<:NTuple{N,Any} where N = LiftingLine{N,T}(models)

"""
    LiftingLine{N}(model)

Construct a lifting line aerodynamic model using `N` instances of `model`.
"""
function LiftingLine{N}(model) where {N,T}
    models = ntuple(i->model, N)
    return LiftingLine(models)
end

# --- Traits --- #

number_of_states(model::LiftingLine) = sum(number_of_states.(model.models))
number_of_inputs(model::LiftingLine) = sum(number_of_inputs.(model.models))
number_of_parameters(model::LiftingLine) = sum(number_of_parameters.(model.models))
inplaceness(::Type{LiftingLine{N,T}}) where {N,T} = InPlace()

function mass_matrix_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(isempty.(mass_matrix_type.(model_types)))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types)))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types)))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(isempty.(state_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types)))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

function input_jacobian_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(isempty.(input_jacobian_type.(model_types)))
        return Empty()
    elseif all(iszero.(input_jacobian_type.(model_types)))
        return Zeros()
    elseif all(isidentity.(input_jacobian_type.(model_types)))
        return Identity()
    elseif all(isconstant.(input_jacobian_type.(model_types)))
        return Constant()
    elseif all(islinear.(input_jacobian_type.(model_types)))
        return Linear()
    else
        return Nonlinear()
    end
end

# --- Methods --- #

function get_rates!(du, model::LiftingLine, u, y, p, t)

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    dus = view.(Ref(du), iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_rates!.(dus, models, us, ys, ps, t)

    return du
end

function get_mass_matrix!(M, model::LiftingLine)

    M .= 0

    models = model.models

    iu = state_indices(model.models)

    Ms = view.(Ref(M), iu, iu)

    get_mass_matrix!.(Ms, models)

    return M
end

function get_mass_matrix!(M, model::LiftingLine, u, y, p, t)

    M .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Ms = view.(Ref(M), iu, iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_mass_matrix!.(Ms, models, us, ys, ps, t)

    return M
end

# --- Performance Overloads --- #

function get_state_jacobian!(J, model::LiftingLine)

    J .= 0

    models = model.models

    iu = state_indices(model.models)

    Js = view.(Ref(J), iu, iu)

    get_state_jacobian!.(Js, models)

    return J
end

function get_state_jacobian!(J, model::LiftingLine, u, y, p, t)

    J .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    Js = view.(Ref(J), iu, iu)
    us = view.(Ref(u), iu)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_state_jacobian!.(Js, models, us, ys, ps, t)

    return J
end

function get_input_jacobian(model::LiftingLine)

    f! = (y, x) -> input_jacobian_product!(y, x, model)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

function get_input_jacobian(model::LiftingLine, λ, d, p, t)

    f! = (y, x) -> input_jacobian_product!(y, x, model, λ, d, p, t)

    M = number_of_states(model)

    N = number_of_inputs(model)

    Jy = LinearMap(f!, M, N; ismutating=true)

    return Jy
end

# --- Unit Testing Methods --- #

# TODO

# function get_lhs(model::LiftingLine, du, u, y, p, t)
#
#     models = model.models
#
#     iu = state_indices(models)
#     iy = input_indices(models)
#     ip = parameter_indices(models)
#
#     dus = view.(Ref(du), iu)
#     us = view.(Ref(u), iu)
#     ys = view.(Ref(y), iy)
#     ps = view.(Ref(p), ip)
#
#     return vcat(get_lhs.(models, dus, us, ys, ps, t)...)
# end

# --- Geometrically Exact Beam Theory Coupling --- #

"""
    couple_models(aero::LiftingLine, stru::GEBT)

Create an aerostructural model using a lifting line aerodynamic model coupled
with a geometrically exact beam theory model.  This model introduces the
freestream velocity components ``\\begin{bmatrix} V_x & V_y & V_z \\end{bmatrix}^T``
as additional parameters.  Also note that when using this model, local reference
frames for each beam should be oriented with the x-axis into the beam, y-axis
forward along the chord-line, and z-axis in the airfoil normal direction.
"""
couple_models(aero::LiftingLine, stru::GEBT)

# --- traits --- #

function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(mass_matrix_type.(model_types, Ref(TypicalSection))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types, Ref(TypicalSection))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types, Ref(TypicalSection))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types, Ref(TypicalSection))))
        return Constant()
    elseif all(islinear.(mass_matrix_type.(model_types, Ref(TypicalSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(state_jacobian_type.(model_types, Ref(TypicalSection))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types, Ref(TypicalSection))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types, Ref(TypicalSection))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types, Ref(TypicalSection))))
        return Constant()
    elseif all(islinear.(state_jacobian_type.(model_types, Ref(TypicalSection))))
        return Linear()
    else
        return Nonlinear()
    end
end

number_of_parameters(::Type{<:LiftingLine}, ::Type{<:GEBT}) = 3

# --- methods --- #

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract freestream velocity (from additional parameters)
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])

    # extract model constants
    system = stru.system
    assembly = stru.assembly

    # loop through each beam element
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract element aerodynamic states
        λi = SVector{Nλi[i]}(λs[i])

        # get structural state variables corresponding to this element
        states = extract_element_state(system, i, q)

        # transformation matrix from local element frame to global frame

        # NOTE: We assume that the local element frame is oriented with the
        # y-axis in the negative chordwise direction and the z-axis in the
        # (airfoil) normal direction

        Ct = GXBeam.get_C(states.theta)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local element frame to local wind frame

        # NOTE: We assume the local wind frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from global frame to the wind frame
        vi = R*CtCab'*Vinf

        # add velocity due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, states.P, states.H)

        # calculate surface angular velocities in the wind frame
        ωi = R*GXBeam.element_angular_velocity(element, states.P, states.H)

        # Note that vi is the **local freestream velocity** in the wind frame,
        # wheras ωi is the **element angular velocity** in the wind frame.

        # calculate lifting line section inputs
        ui = vcat(λi, vi, ωi) # section state variables
        pi = pas[i] # section parameters
        yi = get_inputs(aero.models[i], LiftingLineSection(), ui, pi, t)

        # extract and save aerodynamic model inputs (in wind frame)
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # extract forces and moments per unit length (in wind frame)
        fi = SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        mi = SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # save force and moment per unit length (in global frame)
        ys[6*(i-1)+1 : 6*(i-1)+3] = CtCab*R'*fi
        ys[6*(i-1)+4 : 6*(i-1)+6] = CtCab*R'*mi

    end

    return y
end

function get_input_mass_matrix!(My, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # start with zero valued mass matrix
    My .= 0

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state variables, inputs, and parameters
    λ = view(u, iλ) # aerodynamic state variables
    pa = view(p, ipa) # aerodynamic parameters
    q = view(u, iq) # structural state variables
    ps = view(p, ips) # structural parameters
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract freestream velocity (from additional parameters)
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])

    # extract model constants
    system = stru.system
    assembly = stru.assembly

    # loop through each lifting line / beam element
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract aerodynamic state variables corresponding to this element
        λi = SVector{Nλi[i]}(λs[i])

        # get structural state variables corresponding to this element
        states = extract_element_state(system, i, q)

        # transformation matrix from local element frame to global frame

        # NOTE: We assume that the local element frame is oriented with the
        # y-axis in the negative chordwise direction and the z-axis in the
        # (airfoil) normal direction

        Ct = GXBeam.get_C(states.theta)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local element frame to local wind frame

        # NOTE: We assume the local wind frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from global frame to the wind frame
        vi = R*CtCab'*Vinf

        # add velocity due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, states.P, states.H)

        dv_dP = -R * element.minv11 * system.mass_scaling
        dv_dH = -R * element.minv12 * system.mass_scaling

        # calculate surface angular velocities in the wind frame
        ωi = R*GXBeam.element_angular_velocity(element, states.P, states.H)

        dω_dP = R * element.minv12' * system.mass_scaling
        dω_dH = R * element.minv22 * system.mass_scaling

        # Note that vi is the **local freestream velocity** in the wind frame,
        # wheras ωi is the **element angular velocity** in the wind frame.

        # calculate lifting line section input mass matrix
        ui = vcat(λi, vi, ωi) # section state variables
        pi = pas[i] # section parameters
        Myi = get_input_mass_matrix(aero.models[i], LiftingLineSection(), ui, pi, t)

        # separate into component mass matrices
        d_dλ = SMatrix{Nyai[i], Nλi[i]}(view(Myi, 1:Nyai[i], 1:Nλi[i]))

        d_dv = SMatrix{Nyai[i], 3}(view(Myi, 1:Nyai[i], Nλi[i]+1 : Nλi[i]+3))
        d_dω = SMatrix{Nyai[i], 3}(view(Myi, 1:Nyai[i], Nλi[i]+4 : Nλi[i]+6))
        d_dP = d_dv * dv_dP + d_dω * dω_dP
        d_dH = d_dv * dv_dH + d_dω * dω_dH

        f_dλ = SMatrix{3, Nλi[i]}(view(Myi, Nyai[i]+1 : Nyai[i]+3, 1:Nλi[i]))
        m_dλ = SMatrix{3, Nλi[i]}(view(Myi, Nyai[i]+4 : Nyai[i]+6, 1:Nλi[i]))

        f_dv = SMatrix{3, 3}(view(Myi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+1 : Nλi[i]+3))
        f_dω = SMatrix{3, 3}(view(Myi, Nyai[i]+1 : Nyai[i]+3, Nλi[i]+4 : Nλi[i]+6))
        f_dP = f_dv * dv_dP + f_dω * dω_dP
        f_dH = f_dv * dv_dH + f_dω * dω_dH

        m_dv = SMatrix{3, 3}(view(Myi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+1 : Nλi[i]+3))
        m_dω = SMatrix{3, 3}(view(Myi, Nyai[i]+4 : Nyai[i]+6, Nλi[i]+4 : Nλi[i]+6))
        m_dP = m_dv * dv_dP + m_dω * dω_dP
        m_dH = m_dv * dv_dH + m_dω * dω_dH

        # save aerodynamic input mass matrix entries (in wind frame)
        icol = system.icol_beam[i]
        My[iyas[i], iλs[i]] = d_dλ
        My[iyas[i], nλ+icol+12:nλ+icol+14] = d_dP
        My[iyas[i], nλ+icol+15:nλ+icol+17] = d_dH

        # save load mass matrix entries (in global frame)
        offset = nya + 6*(i-1)
        My[offset+1 : offset+3, iλs[i]] = CtCab*R'*f_dλ
        My[offset+4 : offset+6, iλs[i]] = CtCab*R'*m_dλ
        My[offset+1 : offset+3, nλ+icol+12:nλ+icol+14] = CtCab*R'*f_dP
        My[offset+1 : offset+3, nλ+icol+15:nλ+icol+17] = CtCab*R'*f_dH
        My[offset+4 : offset+6, nλ+icol+12:nλ+icol+14] = CtCab*R'*m_dP
        My[offset+4 : offset+6, nλ+icol+15:nλ+icol+17] = CtCab*R'*m_dH

    end

    return My
end

# --- performance overloads --- #

# TODO

# --- unit testing methods --- #

function get_inputs_from_state_rates(aero::LiftingLine{N,T}, stru::GEBT, du, u, p, t) where {N,T}

    # initialize input vector
    models = (aero, stru)
    TF = promote_type(eltype(du), eltype(u), eltype(p), typeof(t))
    Nu = number_of_inputs(models)
    y = zeros(TF, Nu)

    # extract number of state variables, inputs, and parameters
    nλ = number_of_states(aero) # number of aerodynamic states
    nq = number_of_states(stru) # number of structural states
    nya = number_of_inputs(aero) # number of aerodynamic inputs
    nys = number_of_inputs(stru) # number of structural inputs
    npa = number_of_parameters(aero) # number of aerodynamic parameters
    nps = number_of_parameters(stru) # number of structural parameters
    Nλi = number_of_states.(T.parameters) # number of aerodynamic states for each section
    Nyai = number_of_inputs.(T.parameters) # number of aerodynamic inputs for each section
    Npai = number_of_parameters.(T.parameters) # number of aerodynamic inputs for each section

    # get indices for state variables, inputs, and parameters
    iλ = 1:nλ # indices of aerodynamic states
    iya = 1:nya # indices of aerodynamic inputs
    ipa = 1:npa # indices of aerodynamic parameters
    iq = nλ + 1 : nλ + nq # indices of structural states
    iys = nya + 1 : nya + nys # indices of structural inputs
    ips = npa + 1 : npa + nps # indices of structural parameters
    iλs = state_indices(aero.models) # indices of aerodynamic states for each section
    iyas = input_indices(aero.models) # indices of aerodynamic inputs for each section
    ipas = parameter_indices(aero.models) # indices of aerodynamic inputs for each section

    # separate state rates, states, inputs, and parameters
    dλ = view(du, iλ) # aerodynamic states rates
    λ = view(u, iλ) # aerodynamic state variables
    ya = view(y, iya) # aerodynamic inputs
    pa = view(p, ipa) # aerodynamic parameters
    dq = view(du, iq) # structural state rates
    q = view(u, iq) # structural state variables
    ys = view(y, iys) # structural inputs
    ps = view(p, ips) # structural parameters
    dλs = view.(Ref(dλ), iλs) # aerodynamic state rates for each section
    λs = view.(Ref(λ), iλs) # aerodynamic state variables for each section
    yas = view.(Ref(ya), iyas) # aerodynamic inputs for each section
    pas = view.(Ref(pa), ipas) # aerodynamic parameters for each section

    # extract freestream velocity (from additional parameters)
    Vinf = SVector(p[npa+nps+1], p[npa+nps+2], p[npa+nps+3])

    # extract model constants
    system = stru.system
    assembly = stru.assembly

    # loop through each lifting line / beam element
    for i = 1:N

        # get current structural element
        element = assembly.elements[i]

        # extract element aerodynamic state rates and states
        dλi = SVector{Nλi[i]}(dλs[i])
        λi = SVector{Nλi[i]}(λs[i])

        # get element structural state rates and states
        states = extract_element_state(system, i, q)
        dstates = extract_element_state(system, i, dq)

        # transformation matrix from local element frame to global frame

        # NOTE: We assume that the local element frame is oriented with the
        # y-axis in the negative chordwise direction and the z-axis in the
        # (airfoil) normal direction

        Ct = GXBeam.get_C(states.theta)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # transformation matrix from local element frame to local wind frame

        # NOTE: We assume the local wind frame is oriented with the x-axis
        # in the chordwise direction and the z-axis in the (airfoil) normal
        # direction.

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]

        # transform freestream velocity from global frame to the wind frame
        vi = R*CtCab'*Vinf

        # add velocity due to surface motion
        vi -= R*GXBeam.element_linear_velocity(element, states.P, states.H)

        # calculate freestream acceleration due to surface motion
        dvi = -R*GXBeam.element_linear_velocity(element, dstates.P, dstates.H)

        # calculate surface angular velocities in the wind frame
        ωi = R*GXBeam.element_angular_velocity(element, states.P, states.H)

        # calculate surface angular accelerations in the wind frame
        dωi = R*GXBeam.element_angular_velocity(element, dstates.P, dstates.H)

        # Note that vi is the **local freestream velocity** in the wind frame,
        # wheras ωi is the **element angular velocity** in the wind frame.

        # calculate lifting line section inputs from state rate contributions
        dui = vcat(dλi, dvi, dωi) # section state rates
        ui = vcat(λi, vi, ωi) # section states
        pi = pas[i] # section parameters
        yi = get_inputs_from_state_rates(aero.models[i], LiftingLineSection(),
            dui, ui, pi, t)

        # extract and save aerodynamic model inputs (in wind frame)
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # extract forces and moments per unit length (in wind frame)
        fi = SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        mi = SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # save force and moment per unit length (in global frame)
        ys[6*(i-1)+1 : 6*(i-1)+3] = CtCab*R'*fi
        ys[6*(i-1)+4 : 6*(i-1)+6] = CtCab*R'*mi
    end

    return y
end

# --- Internal Methods --- #

function input_jacobian_product!(y, x, model::LiftingLine)

    models = model.models

    iy = state_indices(models)
    ix = input_indices(models)

    Jyi = get_input_jacobian.(models)

    yi = view.(Ref(y), iy)
    xi = view.(Ref(x), ix)

    mul!.(yi, Jyi, xi)

    return y
end

function input_jacobian_product!(y, x, model::LiftingLine, λ, d, p, t)

    models = model.models

    iu = state_indices(models)
    id = input_indices(models)
    ip = parameter_indices(models)

    xi = view.(Ref(x), id)
    yi = view.(Ref(y), iu)

    λi = view.(Ref(λ), iu)
    di = view.(Ref(d), id)
    pi = view.(Ref(p), ip)

    Jyi = get_input_jacobian.(models, λi, di, pi, t)

    mul!.(yi, Jyi, xi)

    return y
end
