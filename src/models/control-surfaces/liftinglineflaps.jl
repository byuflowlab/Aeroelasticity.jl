"""
    LiftingLine(submodels)

Construct a lifting line model using the aerodynamic models in `models`. State variables, 
inputs, and parameters correspond to the state variables, inputs, and parameters of each of 
the cross sections concatenated
"""
function LiftingLine(submodels)

    # residual
    f = (resid, dx, x, y, p, t) -> liftingline_residual!(resid, dx, x, y, p, t; submodels = submodels)

    # dimensions
    nx = sum(number_of_states.(submodels))
    ny = sum(number_of_inputs.(submodels))
    np = sum(number_of_parameters.(submodels))

    # rate jacobian
    if all(isempty.(rate_jacobian_type.(submodels)))
        ratejac = Empty()
    elseif all(iszero.(rate_jacobian_type.(submodels)))
        ratejac = Zeros()
    elseif all(isidentity.(rate_jacobian_type.(submodels)))
        ratejac = Identity()
    elseif all(isinvariant.(rate_jacobian_type.(submodels)))
        ratejac = Invariant(liftingline_rate_jacobian(; submodels))
    elseif all(isconstant.(rate_jacobian_type.(submodels)))
        ratejac = Constant((J, p) -> liftingline_rate_jacobian!(J, p; submodels))
    elseif all(islinear.(rate_jacobian_type.(submodels)))
        ratejac = Linear(
            (J, dx, x, y, p, t) -> liftingline_rate_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        ratejac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_rate_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # state jacobian
    if all(isempty.(state_jacobian_type.(models)))
        statejac = Empty()
    elseif all(iszero.(state_jacobian_type.(models)))
        statejac = Zeros()
    elseif all(isidentity.(state_jacobian_type.(models)))
        statejac = Identity()
    elseif all(isinvariant.(state_jacobian_type.(models)))
        statejac = Invariant(liftingline_state_jacobian(; submodels))
    elseif all(isconstant.(state_jacobian_type.(models)))
        statejac = Constant((J, p) -> liftingline_state_jacobian!(J, p; submodels))
    elseif all(islinear.(state_jacobian_type.(submodels)))
        statejac = Linear(
            (J, dx, x, y, p, t) -> liftingline_state_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        statejac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_state_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # input jacobian
    if all(isempty.(input_jacobian_type.(models)))
        inputjac = Empty()
    elseif all(iszero.(input_jacobian_type.(models)))
        inputjac = Zeros()
    elseif all(isidentity.(input_jacobian_type.(models)))
        inputjac = Identity()
    elseif all(isinvariant.(input_jacobian_type.(models)))
        inputjac = Invariant(liftingline_input_jacobian(; submodels))
    elseif all(isconstant.(input_jacobian_type.(models)))
        inputjac = Constant((J, p) -> liftingline_input_jacobian!(J, p; submodels))
    elseif all(islinear.(input_jacobian_type.(submodels)))
        inputjac = Linear(
            (J, dx, x, y, p, t) -> liftingline_input_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        inputjac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_input_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # parameter jacobian
    if all(isempty.(parameter_jacobian_type.(models)))
        paramjac = Empty()
    elseif all(iszero.(parameter_jacobian_type.(models)))
        paramjac = Zeros()
    elseif all(isidentity.(parameter_jacobian_type.(models)))
        paramjac = Identity()
    elseif all(isinvariant.(parameter_jacobian_type.(models)))
        paramjac = Invariant(liftingline_parameter_jacobian(; submodels))
    elseif all(isconstant.(parameter_jacobian_type.(models)))
        paramjac = Constant((J, p) -> liftingline_parameter_jacobian!(J, p; submodels))
    elseif all(islinear.(parameter_jacobian_type.(submodels)))
        paramjac = Linear(
            (J, dx, x, y, p, t) -> liftingline_parameter_jacobian!(J, dx, x, y, p, t; submodels)
        )
    else
        paramjac = Nonlinear(
            (J, dx, x, y, p, t) -> liftingline_parameter_jacobian!(J, dx, x, y, p, t; submodels)
        )
    end

    # time gradient
    if all(isempty.(time_gradient_type.(models)))
        tgrad = Empty()
    elseif all(iszero.(time_gradient_type.(models)))
        tgrad = Zeros()
    elseif all(isinvariant.(time_gradient_type.(models)))
        tgrad = Invariant(liftingline_time_gradient(; submodels))
    elseif all(isconstant.(time_gradient_type.(models)))
        tgrad = Constant((dT, p) -> liftingline_time_gradient!(dT, p; submodels))
    elseif all(islinear.(time_gradient_type.(models)))
        paramjac = Linear(
            (dT, dx, x, y, p, t) -> liftingline_time_gradient!(dT, dx, x, y, p, t; submodels)
        )
    else
        paramjac = Nonlinear(
            (dT, dx, x, y, p, t) -> liftingline_time_gradient!(dT, dx, x, y, p, t; submodels)
        )
    end
    
    # input/output functions
    setstate = (x; kwargs...) -> liftingline_setstate!(x; submodels, kwargs...)
    setinput = (y; kwargs...) -> liftingline_setinput!(y; submodels, kwargs...)
    setparam = (p; kwargs...) -> liftingline_setparam!(p; submodels, kwargs...)
    sepstate = (x) -> liftingline_sepstate(x; submodels)
    sepinput = (y) -> liftingline_sepinput(y; submodels)
    sepparam = (p) -> liftingline_sepparam(p; submodels)

    return Model{true}(fresid, nx, ny, np;
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

function liftingline_rate_jacobian(; submodels)
    Jii = get_rate_jacobian.(submodels)
    return block_diagonal(Jii)
end

function liftingline_rate_jacobian!(J, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)
    
    get_rate_jacobian!.(Js, submodels, ps)
    
    return J
end

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

function liftingline_state_jacobian(; submodels)
    Jii = get_state_jacobian.(submodels)
    return block_diagonal(Jii)
end

function liftingline_state_jacobian!(J, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ix)
    ps = view.(Ref(p), ip)
    
    get_state_jacobian!.(Js, submodels, ps)
    
    return J
end

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

function liftingline_input_jacobian(; submodels)
    Jyii = get_input_jacobian.(submodels)
    return block_diagonal(Jyii)
end

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

function liftingline_parameter_jacobian(; submodels)
    Jii = get_parameter_jacobian.(submodels)
    return block_diagonal(Jii)
end

function liftingline_parameter_jacobian!(J, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    iy = input_indices(submodels)
    ip = parameter_indices(submodels)
    
    Js = view.(Ref(J), ix, ip)
    ps = view.(Ref(p), ip)
    
    get_parameter_jacobian!.(Js, submodels, ps)
    
    return J
end

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

liftingline_time_gradient(; submodels) = vcat(get_time_gradient.(submodels)...)

function liftingline_time_gradient!(dT, p; submodels)

    J .= 0
    
    ix = state_indices(submodels)
    ip = parameter_indices(submodels)
    
    dTs = view.(Ref(dT), ix)
    ps = view.(Ref(p), ip)
    
    get_time_gradient!.(dTs, submodels, ps)
    
    return dT
end

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

function liftingline_setstate!(x; section_states, submodels)
    section_indices = state_indices(submodels)
    vxs = view.(Ref(x), section_indices)
    bfn! = (x, model, states) -> set_states!(x, model; states...)
    bfn!.(vxs, submodels, section_states)    
    return x
end

function liftingline_setinput!(y; section_inputs, submodels)
    section_indices = input_indices(submodels)
    vys = view.(Ref(y), section_indices)
    bfn! = (y, model, inputs) -> set_inputs!(y, model; inputs...)
    bfn!.(vys, submodels, section_inputs)
    return y
end

function liftingline_setparam!(p; section_parameters, submodels)
    section_indices = parameter_indices(submodels)
    vps = view.(Ref(p), section_indices)
    bfn! = (p, model, parameters) -> set_parameters!(p, model; parameters...)
    bfn!.(vps, submodels, section_parameters)
    return p
end

function liftingline_sepstate(x; submodels)
    section_indices = state_indices(submodels)
    vxs = view.(Ref(x), section_indices)
    return (section_states = separate_states.(submodels, vxs),)
end

function liftingline_sepinput(y; submodels) 
    section_indices = input_indices(submodels)
    vys = view.(Ref(y), section_indices)
    return (section_inputs = separate_inputs.(submodels, vys),)
end

function liftingline_sepparam(p; submodels)
    section_indices = parameter_indices(submodels)
    vps = view.(Ref(p), section_indices)
    return (section_parameters = separate_parameters.(submodels, vps),)
end