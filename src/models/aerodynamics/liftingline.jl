"""
    LiftingLine(submodels)

Construct a lifting line model using the aerodynamic models in `models`. State variables, 
inputs, and parameters correspond to the state variables, inputs, and parameters of each of 
the cross sections concatenated
"""
function LiftingLine(submodels)

    # model constants (which may be used when this model is coupled with other models)
    constants = submodels

    # residual function
    fresid = (resid, dx, x, y, p, t) -> liftingline_residual!(resid, dx, x, y, p, t; submodels = submodels)

    # number of state, input, and parameters (use Val(N) to use inferrable dimensions)
    nx = sum(number_of_states.(submodels))
    ny = sum(number_of_inputs.(submodels))
    np = sum(number_of_parameters.(submodels))

    # rate jacobian definition (which depends on submodel rate jacobian definitions)
    rate_jacobians = getproperty.(submodels, :ratejac)
    if all(isempty.(rate_jacobians))
        ratejac = Empty()
    elseif all(iszero.(rate_jacobians))
        ratejac = Zeros()
    elseif all(isidentity.(rate_jacobians))
        ratejac = Identity()
    elseif all(isinvariant.(rate_jacobians))
        ratejac = Invariant(liftingline_rate_jacobian(; submodels))
    elseif all(isconstant.(rate_jacobians))
        ratejac = Constant((J, p) -> liftingline_rate_jacobian!(J, p; submodels))
    elseif all(islinear.(rate_jacobians))
        ratejac = Linear(
            (J, dx, x, y, p, t) -> liftingline_rate_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        ratejac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_rate_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # state jacobian definition (which depends on submodel state jacobian definitions)
    state_jacobians = getproperty.(submodels, :statejac)
    if all(isempty.(state_jacobians))
        statejac = Empty()
    elseif all(iszero.(state_jacobians))
        statejac = Zeros()
    elseif all(isidentity.(state_jacobians))
        statejac = Identity()
    elseif all(isinvariant.(state_jacobians))
        statejac = Invariant(liftingline_state_jacobian(; submodels))
    elseif all(isconstant.(state_jacobians))
        statejac = Constant((J, p) -> liftingline_state_jacobian!(J, p; submodels))
    elseif all(islinear.(state_jacobians))
        statejac = Linear(
            (J, dx, x, y, p, t) -> liftingline_state_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        statejac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_state_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # input jacobian definition (which depends on submodel input jacobian definitions)
    input_jacobians = getproperty.(submodels, :inputjac)
    if all(isempty.(input_jacobians))
        inputjac = Empty()
    elseif all(iszero.(input_jacobians))
        inputjac = Zeros()
    elseif all(isidentity.(input_jacobians))
        inputjac = Identity()
    elseif all(isinvariant.(input_jacobians))
        inputjac = Invariant(liftingline_input_jacobian(; submodels))
    elseif all(isconstant.(input_jacobians))
        inputjac = Constant((J, p) -> liftingline_input_jacobian!(J, p; submodels))
    elseif all(islinear.(input_jacobians))
        inputjac = Linear(
            (J, dx, x, y, p, t) -> liftingline_input_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        inputjac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_input_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # parameter jacobian definition (which depends on submodel parameter jacobian definitions)
    parameter_jacobians = getproperty.(submodels, :paramjac)
    if all(isempty.(parameter_jacobians))
        paramjac = Empty()
    elseif all(iszero.(parameter_jacobians))
        paramjac = Zeros()
    elseif all(isidentity.(parameter_jacobians))
        paramjac = Identity()
    elseif all(isinvariant.(parameter_jacobians))
        paramjac = Invariant(liftingline_parameter_jacobian(; submodels))
    elseif all(isconstant.(parameter_jacobians))
        paramjac = Constant((J, p) -> liftingline_parameter_jacobian!(J, p; submodels))
    elseif all(islinear.(parameter_jacobians))
        paramjac = Linear(
            (J, dx, x, y, p, t) -> liftingline_parameter_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        paramjac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_parameter_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # time gradient definition (which depends on submodel time gradient definitions)
    time_gradients = getproperty.(submodels, :tgrad)
    if all(isempty.(time_gradients))
        tgrad = Empty()
    elseif all(iszero.(time_gradients))
        tgrad = Zeros()
    elseif all(isinvariant.(time_gradients))
        tgrad = Invariant(liftingline_time_gradient(; submodels))
    elseif all(isconstant.(time_gradients))
        tgrad = Constant((dT, p) -> liftingline_time_gradient!(dT, p; submodels))
    elseif all(islinear.(time_gradients))
        paramjac = Linear(
            (dT, dx, x, y, p, t) -> liftingline_time_gradient!(dT, dx, x, y, p, t; submodels)
        )
    else
        paramjac = Nonlinear(
            (dT, dx, x, y, p, t) -> liftingline_time_gradient!(dT, dx, x, y, p, t; submodels)
        )
    end
    
    # convenience functions for setting states, inputs, and parameters
    setstate = (x; kwargs...) -> liftingline_setstate!(x; submodels, kwargs...)
    setinput = (y; kwargs...) -> liftingline_setinput!(y; submodels, kwargs...)
    setparam = (p; kwargs...) -> liftingline_setparam!(p; submodels, kwargs...)
    
    # convenience functions for separating states, inputs, and parameters
    sepstate = (x) -> liftingline_sepstate(x; submodels)
    sepinput = (y) -> liftingline_sepinput(y; submodels)
    sepparam = (p) -> liftingline_sepparam(p; submodels)

    # model definition
    return Model{true}(fresid, nx, ny, np;
        constants = constants,
        ratejac = ratejac,
        statejac = statejac,
        inputjac = inputjac,
        paramjac = paramjac,
        tgrad = tgrad,
        setstate = setstate,
        setinput = setinput,
        setparam = setparam,
        sepstate = sepstate,
        sepinput = sepinput,
        sepparam = sepparam,
    )
end

# --- Internal Methods for this Model --- #

# residual function
function liftingline_residual!(resid, dx, x, y, p, t; submodels)

    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)

    resids = view.(Ref(resid), ix)
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)

    get_residual!.(resids, submodels, dxs, xs, ys, ps, t)

    return resid
end

# rate jacobian function (for an invariant rate jacobian)
function liftingline_rate_jacobian(; submodels)
    
    nx = sum(number_of_states.(submodels))

    J = zeros(nx, nx)
    
    ix = state_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    
    get_rate_jacobian!.(Js, submodels)
    
    return J
end

# rate jacobian function (for a constant rate jacobian)
function liftingline_rate_jacobian!(J, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)
    
    get_rate_jacobian!.(Js, submodels, ps)
    
    return J
end

# rate jacobian function (for a time-dependent rate jacobian)
function liftingline_rate_jacobian!(J, dx, x, y, p, t; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_rate_jacobian!.(Js, submodels, dxs, xs, ys, ps, t)
    
    return J
end

# state jacobian function (for an invariant state jacobian)
function liftingline_state_jacobian(; submodels)

    nx = sum(number_of_states.(submodels))

    J = zeros(nx, nx)
    
    ix = state_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    
    get_state_jacobian!.(Js, submodels)

    return J
end

# state jacobian function (for a constant state jacobian)
function liftingline_state_jacobian!(J, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)
    
    get_state_jacobian!.(Js, submodels, ps)
    
    return J
end

# state jacobian function (for a time-dependent state jacobian)
function liftingline_state_jacobian!(J, dx, x, y, p, t; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_state_jacobian!.(Js, submodels, dxs, xs, ys, ps, t)
    
    return J
end

# input jacobian function (for an invariant input jacobian)
function liftingline_input_jacobian(; submodels)

    nx = sum(number_of_states.(submodels))
    ny = sum(number_of_inputs.(submodels))

    J = zeros(nx, ny)
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    
    Js = view.(Ref(J), ix, iy)
    
    get_input_jacobian!.(Js, submodels)

    return J
end

# input jacobian function (for a constant input jacobian)
function liftingline_input_jacobian!(J, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, iy)
    ps = view.(Ref(p), ip)
    
    get_input_jacobian!.(Js, submodels, ps)
    
    return J
end

# input jacobian function (for a time-dependent input jacobian)
function liftingline_input_jacobian!(J, dx, x, y, p, t; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, iy)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_input_jacobian!.(Js, submodels, dxs, xs, ys, ps, t)
    
    return J
end

# parameter jacobian function (for an invariant parameter jacobian)
function liftingline_parameter_jacobian(; submodels)
    
    nx = sum(number_of_states.(submodels))
    np = sum(number_of_parameters.(submodels))

    J = zeros(nx, np)
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ip)
    
    get_parameter_jacobian!.(Js, submodels)

    return J
end

# parameter jacobian function (for a constant parameter jacobian)
function liftingline_parameter_jacobian!(J, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ip)
    ps = view.(Ref(p), ip)
    
    get_parameter_jacobian!.(Js, submodels, ps)
    
    return J
end

# parameter jacobian function (for a time-dependent parameter jacobian)
function liftingline_parameter_jacobian!(J, dx, x, y, p, t; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ip)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_parameter_jacobian!.(Js, submodels, dxs, xs, ys, ps, t)
    
    return J
end

# time gradient function (for an invariant time gradient)
liftingline_time_gradient(; submodels) = vcat(get_time_gradient.(submodels)...)

# time gradient function (for a constant time gradient)
function liftingline_time_gradient!(dT, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    dTs = view.(Ref(dT), ix)
    ps = view.(Ref(p), ip)
    
    get_time_gradient!.(dTs, submodels, ps)
    
    return dT
end

# time gradient function (for a time-dependent time gradient)
function liftingline_time_gradient!(dT, dx, x, y, p, t; submodels)

    dT .= 0
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)
    
    dTs = view.(Ref(dT), ix)
    
    dxs = view.(Ref(dx), ix)
    xs = view.(Ref(x), ix)
    ys = view.(Ref(y), iy)
    ps = view.(Ref(p), ip)
    
    get_time_gradient!.(dTs, submodels, dxs, xs, ys, ps, t)
    
    return J
end

# convenience function for defining this model's state vector
function liftingline_setstate!(x; sections, submodels)
    section_indices = state_indices(submodels)
    vxs = view.(Ref(x), section_indices)
    bfn! = (x, model, states) -> set_states!(x, model; states...)
    bfn!.(vxs, submodels, sections)    
    return x
end

# convenience function for defining this model's input vector
function liftingline_setinput!(y; sections, submodels)
    section_indices = input_indices(submodels)
    vys = view.(Ref(y), section_indices)
    bfn! = (y, model, inputs) -> set_inputs!(y, model; inputs...)
    bfn!.(vys, submodels, sections)
    return y
end

# convenience function for defining this model's parameter vector
function liftingline_setparam!(p; sections, submodels)
    section_indices = parameter_indices(submodels)
    vps = view.(Ref(p), section_indices)
    bfn! = (p, model, parameters) -> set_parameters!(p, model; parameters...)
    bfn!.(vps, submodels, sections)
    return p
end

# convenience function for separating this model's state vector
function liftingline_sepstate(x; submodels)
    section_indices = state_indices(submodels)
    vxs = view.(Ref(x), section_indices)
    return (sections = separate_states.(submodels, vxs),)
end

# convenience function for separating this model's input vector
function liftingline_sepinput(y; submodels) 
    section_indices = input_indices(submodels)
    vys = view.(Ref(y), section_indices)
    return (sections = separate_inputs.(submodels, vys),)
end

# convenience function for separating this model's parameter vector
function liftingline_sepparam(p; submodels)
    section_indices = parameter_indices(submodels)
    vps = view.(Ref(p), section_indices)
    return (sections = separate_parameters.(submodels, vps),)
end

# --- Internal Methods for Couplings with this Model --- #

# local section velocities/accelerations
function liftingline_section_velocities(x, e1, e2, e3, V, Ω, a, α)

    # transformation from body to local frame
    R = [e1 e2 e3]'

    # linear and angular freestream velocity
    vi = -R*(V + cross(Ω, x))
    ωi = R*Ω

    # linear and angular freestream acceleration
    ai = -R*(a + cross(α, x))
    αi = R*α

    return vi, ωi, ai, αi
end