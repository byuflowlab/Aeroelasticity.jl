"""
    LiftingLine{N,T} <: AbstractModel

Lifting line model with `N` cross sections, using the aerodynamic models in `T`.
State variables, inputs, and parameters correspond to the state variables, inputs,
and parameters of each of the cross sections concatenated.
"""
struct LiftingLine{N,T} <: AbstractModel
    models::T
end

# --- Constructors --- #

"""
    LiftingLine(models)

Construct a lifting line aerodynamic model given a tuple of aerodynamic models.
"""
LiftingLine(models::NTuple{N,T}) where {N,T} = LiftingLine{N,T}(models)

LiftingLine{N}(models::NTuple{N,T}) where {N,T} = LiftingLine{N,T}(models)

"""
    LiftingLine{N}(model)

Construct a lifting line aerodynamic model using `N` instances of `model`.
"""
function LiftingLine{N}(model) where {N,T}
    models = ntuple(i->model, N)
    return LiftingLine(models)
end

# --- Traits --- #
number_of_states(model::LiftingLine) = number_of_states(model.models)
number_of_inputs(model::LiftingLine) = number_of_inputs(model.models)
number_of_parameters(model::LiftingLine) = number_of_parameters(model.models)
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
    else
        return Varying()
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
    else
        return Varying()
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
    else
        return Varying()
    end
end

function input_dependence_type(::Type{LiftingLine{N,T}}) where {N,T}
    model_types = (T.parameters...,)
    if all(_linear_input_dependence.(model_types))
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

    get_mass_matrix!.(Ms, models, us, ys, ps)

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

    get_rates!.(dus, models, us, ys, ps)

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

    get_state_jacobian!.(Js, models, us, ys, ps)

    return J
end

function get_input_jacobian!(Jy, model::LiftingLine)

    models = model.models

    iu = state_indices(model.models)
    iy = input_indices(model.models)

    Jys = view.(Ref(Jy), iu, iy)

    get_state_jacobian!.(Jys, models)

    return Jy
end

function get_input_jacobian!(Jy, model::LiftingLine, u, y, p, t)

    Jy .= 0

    models = model.models

    iu = state_indices(models)
    iy = input_indices(models)
    ip = parameter_indices(models)

    us = getindex.(Ref(u), iu)
    ys = getindex.(Ref(y), iy)
    ps = getindex.(Ref(p), ip)

    Jys = view.(Ref(Jy), iu, iy)

    get_state_jacobian!.(Jys, models, us, ys, ps)

    return Jy
end

# TODO: Add parameter jacobian

# --- Coupled Model Properties --- #

# traits
function inplaceness(::Type{LiftingLine{N,T}}, ::Type{<:GEBT}) where {N,T}
    return InPlace()
end

function mass_matrix_type(::Type{LiftingLine{N,T}}, ::Type{GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(mass_matrix_type.(model_types, Ref(GEBT))))
        return Empty()
    elseif all(iszero.(mass_matrix_type.(model_types, Ref(GEBT))))
        return Zeros()
    elseif all(isidentity.(mass_matrix_type.(model_types, Ref(GEBT))))
        return Identity()
    elseif all(isconstant.(mass_matrix_type.(model_types, Ref(GEBT))))
        return Constant()
    else
        return Varying()
    end
end

function state_jacobian_type(::Type{LiftingLine{N,T}}, ::Type{GEBT}) where {N,T}
    model_types = T.parameters
    if all(isempty.(state_jacobian_type.(model_types, Ref(GEBT))))
        return Empty()
    elseif all(iszero.(state_jacobian_type.(model_types, Ref(GEBT))))
        return Zeros()
    elseif all(isidentity.(state_jacobian_type.(model_types, Ref(GEBT))))
        return Identity()
    elseif all(isconstant.(state_jacobian_type.(model_types, Ref(GEBT))))
        return Constant()
    else
        return Varying()
    end
end

function get_input_mass_matrix!(My, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # separate aerodynamic and structural state variables
    iλ = state_indices(aero)
    iq = state_indices(stru)
    λ = view(u, iλ)
    q = view(u, iq)
    λs = view.(Ref(λ), state_indices(aero.models))

    # separate aerodynamic and structural parameters
    ipa = parameter_indices(aero)
    ips = parameter_indices(stru)
    pa = view(p, ip_a)
    ps = view(p, ip_s)
    pas = view.(Ref(pa), parameter_indices(aero.models))

    # separate aerodynamic and structural state rates variables
    dλ = view(du, iλ)
    dq = view(du, iq)
    dλs = view.(Ref(dλ), state_indices(aero.models))

    # extract model properties
    system = stru.system
    assembly = stru.assembly

    # loop through each lifting line / beam element
    for i = 1:N
        # aerodynamic model
        aerodynamic_model = aero.models[i]

        # beam element properties and state variables
        element = assembly.elements[i]
        states = extract_element_state(system, i; x = q)

        # beam element deflections and velocities
        dP_dPi = Diagonal((@SVector ones(3)))
        dH_dHi = Diagonal((@SVector ones(3))) .* system.mass_scaling

        dv_dPi = element.minv11 .* system.mass_scaling
        dv_dHi = element.minv12 .* system.mass_scaling
        dω_dPi = element.minv21 .* system.mass_scaling
        dω_dHi = element.minv22 .* system.mass_scaling

        # structural state variables
        dhdot = dv_dPi[3] # plunge acceleration
        dθdot = dω_dHi[2] # pitch acceleration
        qi = SVector(h, θ, hdot, θdot)

        # combined state variables
        ui = vcat(λi, qi)

        # aerodynamic parameters
        pai = pas[i]

        # structural parameters (currently unused by all aerodynamic models)


        # coupled parameters
        pi = pai

        # NOTE: We currently only include aerodynamic parameters in the section
        # parameter vector.  This would have to be changed if a 2D model uses
        # structural parameters to calculate coupled inputs

        # combined section inputs
        Myi = get_input_mass_matrix(aero.models[i], TypicalSection(), ui, pi, t)

        # set section aerodynamic inputs
        yas[i] .= view(yi, 1:length(yi)-2)

        # set section structural inputs
        L = yi[end-1]
        M = yi[end]

        # convert lift and moment to distributed loads
        ys[3*(i-1)+1] = 0
        ys[3*(i-1)+2] = 0
        ys[3*(i-1)+3] = L
        ys[3*(i-1)+4] = 0
        ys[3*(i-1)+5] = M + (b/2+a*b)*L
        ys[3*(i-1)+6] = 0
    end

    return y
end

function get_inputs!(y, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # separate aerodynamic and structural state variables
    iλ = state_indices(aero)
    iq = state_indices(stru)
    λ = view(u, iλ)
    q = view(u, iq)
    λs = view.(Ref(λ), state_indices(aero.models))

    # separate aerodynamic and structural inputs
    iya = input_indices(aero)
    iys = input_indices(stru)
    ya = view(y, iy_a)
    ys = view(y, iy_s)
    yas = view.(Ref(ya), input_indices(aero.models))

    # separate aerodynamic and structural parameters
    ipa = parameter_indices(aero)
    ips = parameter_indices(stru)
    pa = view(p, ip_a)
    ps = view(p, ip_s)
    pas = view.(Ref(pa), parameter_indices(aero.models))

    # extract model properties
    system = stru.system
    assembly = stru.assembly

    # loop through each lifting line / beam element
    for i = 1:N
        # aerodynamic model
        aerodynamic_model = aero.models[i]

        # beam element properties and state variables
        element = assembly.elements[i]
        states = extract_element_state(system, i; x = q)

        # beam element deflections
        u_elem = states.u

        # beam element rotation parameters
        θ_elem = states.θ

        # transformation matrix from global to deformed local frame
        Ct = get_C(θ)'
        Cab = beam.Cab
        CtCab = Ct*Cab

        # linear and angular velocities (transform from deformed local to global frame)
        v_elem = CtCab' * GXBeam.element_linear_velocity(element, states.P, states.H)
        ω_elem = CtCab' * GXBeam.element_angular_velocity(element, states.P, states.H)

        # aerodynamic state variables
        λi = λs[i] #TODO: make this statically sized

        # structural state variables
        h = states.u[3] # plunge
        θ = atan(CtCab[3,1], CtCab[1,1]) # pitch (in plane parallel to Y-Z plane)
        hdot = v_elem[3] # plunge rate
        θdot = ω_elem[2] # pitch rate (in plane parallel to Y-Z plane)
        qi = SVector(h, θ, hdot, θdot)

        # combined state variables
        ui = vcat(λi, qi)

        # aerodynamic parameters
        pai = pas[i]

        # structural parameters (currently unused by all aerodynamic models)


        # coupled parameters
        pi = pai

        # NOTE: We currently only include aerodynamic parameters in the section
        # parameter vector.  This would have to be changed if a 2D model uses
        # structural parameters to calculate coupled inputs

        # combined section inputs
        yi = get_inputs(aero.models[i], TypicalSection(), ui, pi, t)

        # set section aerodynamic inputs
        yas[i] .= view(yi, 1:length(yi)-2)

        # set section structural inputs
        L = yi[end-1]
        M = yi[end]

        # convert lift and moment to distributed loads
        ys[3*(i-1)+1] = 0
        ys[3*(i-1)+2] = 0
        ys[3*(i-1)+3] = L
        ys[3*(i-1)+4] = 0
        ys[3*(i-1)+5] = M + (b/2+a*b)*L
        ys[3*(i-1)+6] = 0
    end

    return y
end

function get_input_state_jacobian!(My, aero::LiftingLine{N,T}, stru::GEBT, u, p, t) where {N,T}

    # separate aerodynamic and structural state variables
    iλ = state_indices(aero)
    iq = state_indices(stru)
    λ = view(u, iλ)
    q = view(u, iq)
    λs = view.(Ref(λ), state_indices(aero.models))

    # separate aerodynamic and structural parameters
    ipa = parameter_indices(aero)
    ips = parameter_indices(stru)
    pa = view(p, ip_a)
    ps = view(p, ip_s)
    pas = view.(Ref(pa), parameter_indices(aero.models))

    # separate input mass matrix
    iya = input_indices(aero)
    iys = input_indices(stru)
    ya = view(y, iy_a)
    ys = view(y, iy_s)
    yas = view.(Ref(ya), input_indices(aero.models))

    # extract model properties
    system = stru.system
    assembly = stru.assembly

    # loop through each lifting line / beam element
    for i = 1:N
        # aerodynamic model
        aerodynamic_model = aero.models[i]

        # beam element deflections and velocities
        icol = system.icol_beam[i]
        u_ui = Diagonal((@SVector ones(3)))
        θ_θi = Diagonal((@SVector ones(3)))
        P_Pi = Diagonal((@SVector ones(3))) .* system.mass_scaling
        H_Hi = Diagonal((@SVector ones(3))) .* system.mass_scaling

        v_Pi = element.minv11*P_Pi
        v_Hi = element.minv12*H_Hi
        ω_Pi = element.minv21*P_Pi
        ω_Hi = element.minv22*H_Hi

        # convert rotation parameter to Wiener-Milenkovic parameters
        scaling = rotation_parameter_scaling(theta)
        θ_θi *= scaling

        # aerodynamic state variables
        λi = λs[i] #TODO: make this statically sized

        # structural state variables
        h = states.u[3] # plunge
        θ = # pitch
        hdot = udot_elem[3] # plunge rate
        θdot =  # pitch rate
        qi = SVector(h, θ, hdot, θdot)

        # combined state variables
        ui = vcat(λi, qi)

        # aerodynamic parameters
        pai = pas[i]

        # structural parameters (currently unused by all aerodynamic models)


        # coupled parameters
        pi = pai

        # NOTE: We currently only include aerodynamic parameters in the section
        # parameter vector.  This would have to be changed if a 2D model uses
        # structural parameters to calculate coupled inputs

        # combined section inputs
        Myi = get_input_state_jacobian(aero.models[i], TypicalSection(), ui, pi, t)

        # set section aerodynamic inputs
        My[]
        yas[iys[i], ] = Myi[1:length(yi)-2, :] * qi_

        # set section structural inputs
        L = yi[end-1]
        M = yi[end]

        # convert lift and moment to distributed loads
        ys[3*(i-1)+1] = 0
        ys[3*(i-1)+2] = 0
        ys[3*(i-1)+3] = L
        ys[3*(i-1)+4] = 0
        ys[3*(i-1)+5] = M + (b/2+a*b)*L
        ys[3*(i-1)+6] = 0
    end

    return y
end
