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
``r = \\begin{bmatrix} L' & Y' & D' & Mx, My, Mz \\end{bmatrix}^T``, and zero
parameters.  Two-dimensional aerodynamic models may be extended to
three-dimensional models by coupling with this model.  Note that
this model has no rate equations of its own.
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

# --- Coupled Model Properties --- #

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

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{GEBT}) where {N,T}
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

    iv1 = SVector(1, 2, 3)
    iv2 = SVector(4, 5, 6)

    # aerodynamic state variables, inputs, and parameters
    iλ = state_indices(aero)
    iya = input_indices(aero)
    ipa = parameter_indices(aero)

    λ = view(u, iλ)
    ya = view(y, iy_a)
    pa = view(p, ip_a)

    # structural state variables, inputs, and parameters
    iq = state_indices(stru)
    iys = input_indices(stru)
    ips = parameter_indices(stru)

    q = view(u, iq)
    ys = view(y, iy_s)
    ps = view(p, ip_s)

    # section aerodynamic state variables, inputs, and parameters
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

        # extract aerodynamic model and state variables
        aerodynamic_model = aero.models[i]
        λi = λs[i]

        # extract structural element properties and state variables
        element = assembly.elements[i]
        states = extract_element_state(system, i; x = q)

        # calculate transformation matrix from local to global frame
        Ct = GXBeam.get_C(states.theta)'
        Cab = beam.Cab
        CtCab = Ct*Cab

        # calculate linear and angular velocities
        v_elem = GXBeam.element_linear_velocity(element, states.P, states.H)
        ω_elem = GXBeam.element_angular_velocity(element, states.P, states.H)
        qi = vcat(v_elem, ω_elem)

        dv_dP = element.minv11 * system.mass_scaling
        dv_dH = element.minv12 * system.mass_scaling
        dω_dP = element.minv21 * system.mass_scaling
        dω_dP = element.minv22 * system.mass_scaling

        # calculate section aerodynamic input mass matrices
        d_dλ, d_dqi = section_aerodynamic_input_mass_matrices(aerodynamic_model,
            λi, qi, pas[i], t)

        d_dv = d_dqi[:, iv1]
        d_dω = d_dqi[:, iv2]

        d_dP = d_dv * dv_dP + d_dω * dω_dP
        d_dH = d_dv * dv_dH + d_dω * dω_dH

        # calculate section structural input mass matrices
        r_dλ, r_dqi = section_aerodynamic_load_mass_matrices(aerodynamic_model,
            λi, qi, pas[i], t)

        r_dv = r_dqi[:, iv1]
        r_dω = r_dqi[:, iv2]

        r_dP = r_dv * dv_dP + r_dω * dω_dP
        r_dH = r_dv * dv_dH + r_dω * dω_dH

        # extract force and moment per unit length
        f_dλi = r_dλ[iv1, :]
        f_dP = r_dP[iv1, :]
        f_dH = r_dH[iv1, :]

        m_dλi = r_dλ[iv2, :]
        m_dP = r_dP[iv2, :]
        m_dH = r_dH[iv2, :]

        iP = iq[SVector(icol+12, icol+13, icol+14)]
        iH = iq[SVector(icol+15, icol+16, icol+17)]

        ifi = 6*(i-1) .+ iv1
        imi = 6*(i-1) .+ iv2
        # save results
        My[iyas[i], iλs[i]] .= d_dλi
        My[iyas[i], iP] .= d_dPi
        My[iyas[i], iH] .= d_dHi
        My[ifi, iλs[i]] .= CtCab*f_dλi
        My[imi, iλs[i]] .= CtCab*m_dλi
        My[ifi, iP] .= CtCab*f_dPi
        My[imi, iP] .= CtCab*m_dPi
        My[ifi, iH] .= CtCab*f_dHi
        My[imi, iH] .= CtCab*m_dHi
    end

    return y
end

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # aerodynamic state variables, inputs, and parameters
    iλ = state_indices(aero)
    iya = input_indices(aero)
    ipa = parameter_indices(aero)

    λ = view(u, iλ)
    ya = view(y, iy_a)
    pa = view(p, ip_a)

    # structural state variables, inputs, and parameters
    iq = state_indices(stru)
    iys = input_indices(stru)
    ips = parameter_indices(stru)

    q = view(u, iq)
    ys = view(y, iy_s)
    ps = view(p, ip_s)

    # section aerodynamic state variables, inputs, and parameters
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

        # extract aerodynamic model and state variables
        aerodynamic_model = aero.models[i]
        λi = λs[i]

        # extract structural element properties and state variables
        element = assembly.elements[i]
        states = extract_element_state(system, i; x = q)

        # calculate transformation matrix from local to global frame
        Ct = GXBeam.get_C(states.theta)'
        Cab = beam.Cab
        CtCab = Ct*Cab

        # calculate linear and angular velocity
        v_elem = GXBeam.element_linear_velocity(element, states.P, states.H)
        ω_elem = GXBeam.element_angular_velocity(element, states.P, states.H)
        qi = vcat(v_elem, ω_elem)

        # calculate aerodynamic inputs (in local frame)
        d = section_aerodynamic_inputs(aerodynamic_model, λi, qi, pas[i], t)

        # calculate aerodynamic loads (in local frame)
        r = section_aerodynamic_loads(aerodynamic_model, λi, qi, pas[i], t)

        # extract force and moment per unit length
        f = SVector(r[1], r[2], r[3])
        m = SVector(r[4], r[5], r[6])

        # save results
        ys[iyas[i]] = d # aerodynamic inputs (in local frame)
        ys[6*(i-1)+1 : 6*(i-1)+3] = CtCab*f # aerodynamic loads (in global frame)
        ys[6*(i-1)+4 : 6*(i-1)+6] = CtCab*m # aerodynamic loads (in global frame)

    end

    return y
end
