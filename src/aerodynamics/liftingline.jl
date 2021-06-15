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

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    Ms = view.(Ref(M), iu, iu)

    get_mass_matrix!.(Ms, models, us, ys, ps, t)

    return M
end

function get_rates!(du, model::LiftingLine, u, y, p, t)

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    dus = view.(Ref(du), iu)

    get_rates!.(dus, models, us, ys, ps, t)

    return du
end

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

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    Js = view.(Ref(J), iu, iu)

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

# TODO: Add parameter jacobian

# --- Geometrically Exact Beam Theory --- #

# traits
function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{GEBT}) where {N,T}
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

function get_input_mass_matrix!(My, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # aerodynamic state variables, inputs, and parameters
    nλ = number_of_states(aero)
    nya = number_of_inputs(aero)
    npa = number_of_parameters(aero)

    iλ = 1:nλ
    iya = 1:nya
    ipa = 1:npa

    λ = view(u, iλ)
    ya = view(y, iya)
    pa = view(p, ipa)

    # structural state variables, inputs, and parameters
    nq = number_of_states(stru)
    nys = number_of_inputs(stru)
    nps = number_of_parameters(stru)

    iq = nλ + 1 : nλ + nq
    iys = nya + 1 : nya + nys
    ips = npa + 1 : npa + nps

    q = view(u, iq)
    ys = view(y, iys)
    ps = view(p, ips)

    # aerodynamic section state variables, inputs, and parameters
    Nλi = number_of_states.(T.parameters)
    Nyai = number_of_inputs.(T.parameters)
    Npai = number_of_parameters.(T.parameters)

    iλs = state_indices(aero.models)
    iyas = input_indices(aero.models)
    ipas = parameter_indices(aero.models)

    λs = view.(Ref(λ), iλs)
    yas = view.(Ref(ya), iyas)
    pas = view.(Ref(pa), ipas)

    # structural model properties
    system = stru.system
    assembly = stru.assembly

    # loop through each lifting line / beam element
    for i = 1:N

        # extract aerodynamic state variables
        λi = SVector{Nλi}(λs[i])

        # calculate lifting line section state variables

        # element properties and state variables
        element = assembly.elements[i]
        states = extract_element_state(system, i, q)

        # transformation matrix from local aerodynamic frame to global frame
        Ct = GXBeam.get_C(states.theta)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # linear and angular velocity (in local structural frame)
        velem = GXBeam.element_linear_velocity(element, states.P, states.H)
        ωelem = GXBeam.element_angular_velocity(element, states.P, states.H)

        dv_dP = element.minv11 * system.mass_scaling
        dv_dH = element.minv12 * system.mass_scaling
        dω_dP = element.minv21 * system.mass_scaling
        dω_dP = element.minv22 * system.mass_scaling

        # linear and angular velocity (in local aerodynamic frame)
        vi = SVector(-velem[2], velem[1], velem[3])
        ωi = SVector(-ωelem[2], ωelem[1], ωelem[3])

        R = @SMatrix [0 -1 0; 1 0 0; 0 0 1]
        dv_dP = R*dv_dP
        dv_dH = R*dv_dH
        dω_dP = R*dω_dP
        dω_dP = R*dω_dP

        # state variable vector for this section
        ui = vcat(λi, vi, ωi)

        # parameter vector for this section
        pi = pas[i]

        # input mass matrix for this section
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

        # save aerodynamic input mass matrix entries
        My[iyas[i], iλs[i]] .= d_dλ
        My[iyas[i], nλ+icol+12:nλ+icol+14] .= d_P
        My[iyas[i], nλ+icol+15:nλ+icol+17] .= d_H

        # convert to local structural frame of reference
        f_dλ = R'*f_dλ
        m_dλ = R'*m_dλ
        f_dP = R'*f_dP
        f_dH = R'*f_dH
        m_dP = R'*m_dP
        m_dH = R'*m_dH

        # save load mass matrix entries (in global frame of reference)
        offset = nya + 6*(i-1)
        My[offset+1 : offset+3, iλs[i]] .= CtCab*f_dλ
        My[offset+4 : offset+6, iλs[i]] .= CtCab*m_dλ
        My[offset+1 : offset+3, nλ+icol+12:nλ+icol+14] .= CtCab*f_dP
        My[offset+1 : offset+3, nλ+icol+15:nλ+icol+17] .= CtCab*f_dH
        My[offset+4 : offset+6, nλ+icol+12:nλ+icol+14] .= CtCab*m_dP
        My[offset+4 : offset+6, nλ+icol+15:nλ+icol+17] .= CtCab*m_dH

    end

    return y
end

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # aerodynamic state variables, inputs, and parameters
    nλ = number_of_states(aero)
    nya = number_of_inputs(aero)
    npa = number_of_parameters(aero)

    iλ = 1:nλ
    iya = 1:nya
    ipa = 1:npa

    λ = view(u, iλ)
    ya = view(y, iya)
    pa = view(p, ipa)

    # structural state variables, inputs, and parameters
    nq = number_of_states(stru)
    nys = number_of_inputs(stru)
    nps = number_of_parameters(stru)

    iq = nλ + 1 : nλ + nq
    iys = nya + 1 : nya + nys
    ips = npa + 1 : npa + nps

    q = view(u, iq)
    ys = view(y, iys)
    ps = view(p, ips)

    # aerodynamic section state variables, inputs, and parameters
    Nλi = number_of_states.(T.parameters)
    Nyai = number_of_inputs.(T.parameters)
    Npai = number_of_parameters.(T.parameters)

    iλs = state_indices(aero.models)
    iyas = input_indices(aero.models)
    ipas = parameter_indices(aero.models)

    λs = view.(Ref(λ), iλs)
    yas = view.(Ref(ya), iyas)
    pas = view.(Ref(pa), ipas)

    # structural model properties
    system = stru.system
    assembly = stru.assembly

    # loop through each lifting line / beam element
    for i = 1:N

        # extract aerodynamic state variables
        λi = SVector{Nλi[i]}(λs[i])

        # calculate lifting line section state variables

        # element properties and state variables
        element = assembly.elements[i]
        states = extract_element_state(system, i, q)

        # transformations matrix from local aerodynamic frame to global frame
        Ct = GXBeam.get_C(states.theta)'
        Cab = element.Cab
        CtCab = Ct*Cab

        # linear and angular velocity (in local frame)
        velem = GXBeam.element_linear_velocity(element, states.P, states.H)
        ωelem = GXBeam.element_angular_velocity(element, states.P, states.H)

        # linear and angular velocity experienced by this section
        vi = SVector(-velem[2], velem[1], velem[3])
        ωi = SVector(-ωelem[2], ωelem[1], ωelem[3])

        # state variable vector for this section
        ui = vcat(λi, vi, ωi)

        # parameter vector for this section
        pi = pas[i]

        # inputs for this section
        yi = get_inputs(aero.models[i], LiftingLineSection(), ui, pi, t)

        # save aerodynamic inputs
        for iy = 1:Nyai[i]
            yas[i][iy] = yi[iy]
        end

        # get forces/moments per unit length
        fi = SVector(yi[Nyai[i]+1], yi[Nyai[i]+2], yi[Nyai[i]+3])
        mi = SVector(yi[Nyai[i]+4], yi[Nyai[i]+5], yi[Nyai[i]+6])

        # convert to element frame of reference
        f_elem = SVector(fi[2], -fi[1], fi[3])
        m_elem = SVector(mi[2], -mi[1], mi[3])

        # save loads (in global frame of reference)
        ys[6*(i-1)+1 : 6*(i-1)+3] = CtCab*f_elem
        ys[6*(i-1)+4 : 6*(i-1)+6] = CtCab*m_elem

    end

    return y
end
