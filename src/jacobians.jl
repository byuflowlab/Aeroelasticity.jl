# --- Jacobian Types --- #

abstract type AbstractJacobian end

"""
    Empty <: AbstractJacobian

A gradient/jacobian wrapper which indicates that a given gradient/jacobian is empty.
"""
struct Empty <: AbstractJacobian end

"""
    Zeros <: AbstractJacobian

A gradient/jacobian wrapper which indicates that a given gradient/jacobian is zero.
"""
struct Zeros <: AbstractJacobian end

"""
    Identity <: AbstractJacobian

A gradient/jacobian wrapper which indicates that a given gradient/jacobian is the identity 
matrix.
"""
struct Identity <: AbstractJacobian end

"""
    Invariant <: AbstractJacobian

A gradient/jacobian wrapper which indicates that the value of given gradient/jacobian is 
independent of the state rates, states, inputs, parameters, and time.

The value of the invariant jacobian/gradient may be provided upon construction, otherwise
it will be calculated using automatic differentiation.
"""
struct Invariant{T} <: AbstractJacobian; value::T; end

Invariant() = Invariant(nothing)

"""
    Constant <: AbstractJacobian

A gradient/jacobian wrapper which indicates that the value of given gradient/jacobian is 
independent of the state rates, states, inputs, and time.

The jacobian/gradient will be calculated using automatic differentiation.  Alternatively,
a function may be provided which defines the value of the gradient/jacobian as a function 
of the parameters.
"""
struct Constant{T} <: AbstractJacobian; func::T; end

Constant() = Constant(nothing)

"""
    Linear <: AbstractJacobian

A gradient/jacobian wrapper which indicates that the value of given gradient/jacobian is 
linear with respect to the differentiated variable.

The jacobian/gradient will be calculated using automatic differentiation.  Alternatively,
a function may be provided which defines the value of the gradient/jacobian as a function 
of the state rates, states, inputs, parameters, and time. 
"""
struct Linear{T} <: AbstractJacobian; func::T; end

Linear() = Linear(nothing)

"""
    Nonlinear <: AbstractJacobian

A gradient/jacobian wrapper which indicates that the value of given gradient/jacobian is 
nonlinear with respect to the differentiated variable.

The jacobian/gradient will be calculated using automatic differentiation.  Alternatively,
a function may be provided which defines the value of the gradient/jacobian as a function 
of the state rates, states, inputs, parameters, and time. 
"""
struct Nonlinear{T} <: AbstractJacobian; func::T; end

Nonlinear() = Nonlinear(nothing)

# trait dispatch
isempty(::AbstractJacobian) = false
isempty(::Empty) = true

iszero(::AbstractJacobian) = false
iszero(::Zeros) = true

isidentity(::AbstractJacobian) = false
isidentity(::Identity) = true

isinvariant(::AbstractJacobian) = false
isinvariant(::Empty) = true
isinvariant(::Zeros) = true
isinvariant(::Identity) = true
isinvariant(::Invariant) = true

isconstant(::AbstractJacobian) = false
isconstant(::Empty) = true
isconstant(::Zeros) = true
isconstant(::Identity) = true
isconstant(::Invariant) = true
isconstant(::Constant) = true

islinear(::AbstractJacobian) = false
islinear(::Empty) = true
islinear(::Zeros) = true
islinear(::Identity) = true
islinear(::Invariant) = true
islinear(::Constant) = true
islinear(::Linear) = true

# --- Jacobian Update Functions --- #

for iarg = (1, 2, 3, 4, 5)

    if iarg == 1
        fname = :add_rate_jacobian
        adfunc = :autodiff_jacobian_func
        ncol = :nx
    elseif iarg == 2
        fname = :add_state_jacobian
        adfunc = :autodiff_jacobian_func
        ncol = :nx
    elseif iarg == 3
        fname = :add_input_jacobian
        adfunc = :autodiff_jacobian_func
        ncol = :ny
    elseif iarg == 4
        fname = :add_parameter_jacobian
        adfunc = :autodiff_jacobian_func
        ncol = :np
    elseif iarg == 5
        fname = :add_time_gradient
        adfunc = :autodiff_derivative_func
        ncol = 1
    end

    @eval begin
        # by default, return original rate jacobian
        $(fname)(jac, f, iip, nx, ny, np) = jac

        # invariant jacobian
        function $(fname)(::Invariant{Nothing}, f, iip, nx, ny, np)
               
            # set values arbitrarily
            dx = zeros(val(nx))
            x = zeros(val(nx))
            y = zeros(val(ny))
            p = zeros(val(np))
            t = 0.0

            if iip
                # initialize jacobian
                J = zeros(nx, $ncol)
                # get in-place function
                fJ = $(adfunc)(f, val(nx), $iarg)
                # calculate in-place jacobian
                fJ(J, dx, x, y, p, t)
            else # out-of-place
                # get out-of-place function
                fJ = $(adfunc)(f, $iarg)
                # calculate out-of-place jacobian 
                J = fJ(dx, x, y, p, t)
            end

            return Invariant(J)
        end

        # constant jacobian
        function $(fname)(::Constant{Nothing}, f, iip, nx, ny, np)
    
            # set values arbitrarily
            dx = zeros(nx)
            x = zeros(nx)
            y = zeros(ny)
            t = 0.0

            if iip
                # get in-place function
                fJ = $(adfunc)(f, val(nx), $iarg)
                # get in-place jacobian function
                fJ = (J, p) -> fJ(J, dx, x, y, p, t)
            else
                # get out-of-place function
                fJ = $(adfunc)(f, $iarg)
                # get out-of-place jacobian function
                fJ = (p) -> fJ(dx, x, y, p, t)

            end
        
            return Constant(fJ)
        end

        # out-of-place or in-place linear rate jacobian
        function $(fname)(::Linear{Nothing}, f, iip, nx, ny, np)
            if iip
                fJ = $(adfunc)(f, nx, $iarg)
            else
                fJ = $(adfunc)(f, $iarg)
            end
            return Linear(fJ)
        end

        # out-of-place or in-place nonlinear rate jacobian
        function $(fname)(::Nonlinear{Nothing}, f, iip, nx, ny, np)
            if iip
                fJ = $(adfunc)(f, nx, $iarg)
            else
                fJ = $(adfunc)(f, $iarg)
            end
            return Nonlinear(fJ)
        end

    end
end

for iarg = (1, 2, 3, 4)

    if iarg == 1
        fname = :add_coupling_rate_jacobian
        adfunc = :autodiff_jacobian_func
        ncol = :nx
    elseif iarg == 2
        fname = :add_coupling_state_jacobian
        adfunc = :autodiff_jacobian_func
        ncol = :nx
    elseif iarg == 3
        fname = :add_coupling_parameter_jacobian
        adfunc = :autodiff_jacobian_func
        ncol = :np
    elseif iarg == 4
        fname = :add_coupling_time_gradient
        adfunc = :autodiff_derivative_func
        ncol = 1
    end

    @eval begin
        # by default, return original rate jacobian
        $(fname)(jac, f, iip, nx, ny, np) = jac

        # invariant jacobian
        function $(fname)(::Invariant{Nothing}, f, iip, nx, ny, np)
               
            # set values arbitrarily
            dx = zeros(val(nx))
            x = zeros(val(nx))
            p = zeros(val(np))
            t = 0.0

            if iip
                # initialize jacobian
                J = zeros(nx, $ncol)
                # get in-place function
                fJ = $(adfunc)(f, val(ny), $iarg)
                # calculate in-place jacobian
                fJ(J, dx, x, p, t)
            else # out-of-place
                # get out-of-place function
                fJ = $(adfunc)(f, $iarg)
                # calculate out-of-place jacobian 
                J = fJ(dx, x, p, t)
            end

            return Invariant(J)
        end

        # constant jacobian
        function $(fname)(::Constant{Nothing}, f, iip, nx, ny, np)
    
            # set values arbitrarily
            dx = zeros(nx)
            x = zeros(nx)
            t = 0.0

            if iip
                # get in-place function
                fJ = $(adfunc)(f, val(ny), $iarg)
                # get in-place jacobian function
                fJ = (J, p) -> fJ(J, dx, x, p, t)
            else
                # get out-of-place function
                fJ = $(adfunc)(f, $iarg)
                # get out-of-place jacobian function
                fJ = (p) -> fJ(dx, x, p, t)

            end
        
            return Constant(fJ)
        end

        # out-of-place or in-place linear rate jacobian
        function $(fname)(::Linear{Nothing}, f, iip, nx, ny, np)
            if iip
                fJ = $(adfunc)(f, val(ny), $iarg)
            else
                fJ = $(adfunc)(f, $iarg)
            end
            return Linear(fJ)
        end

        # out-of-place or in-place nonlinear rate jacobian
        function $(fname)(::Nonlinear{Nothing}, f, iip, nx, ny, np)
            if iip
                fJ = $(adfunc)(f, val(ny), $iarg)
            else
                fJ = $(adfunc)(f, $iarg)
            end
            return Nonlinear(fJ)
        end

    end
end

# --- Jacobian Evaluation --- #

# out-of-place, statically sized

function jacobian(::Empty, ::Val{NX}, ::Val{NY}, args...) where {NX, NY}
    return SMatrix{NX, NY, Float64}()
end

function jacobian(::Zeros, ::Val{NX}, ::Val{NY}, args...) where {NX, NY}
    return zeros(SMatrix{NX, NY, Float64})
end

function jacobian(::Identity, ::Val{NX}, ::Val{NY}, args...) where {NX, NY}
    return SMatrix{NX, NY}(I)
end

function jacobian(jac::Invariant, ::Val{NX}, ::Val{NY}, args...) where {NX, NY}
    return SMatrix{NX, NY}(jac.value)
end

function jacobian(jac::Constant, ::Val{NX}, ::Val{NY}, p) where {NX, NY}
    return SMatrix{NX, NY}(jac.func(p))
end

function jacobian(jac::Constant, ::Val{NX}, ::Val{NY}, dx, x, p, t) where {NX, NY}
    return SMatrix{NX, NY}(jac.func(p))
end

function jacobian(jac::Constant, ::Val{NX}, ::Val{NY}, dx, x, y, p, t) where {NX, NY}
    return SMatrix{NX, NY}(jac.func(p))
end

function jacobian(jac::Union{Linear, Nonlinear}, ::Val{NX}, ::Val{NY}, args...) where {NX, NY}
    return SMatrix{NX, NY}(jac.func(args...))
end

# out-of-place, dynamically sized

function jacobian(::Empty, nx, ny, args...)
    return Matrix{Float64}(nx, ny)
end

function jacobian(::Zeros, nx, ny, args...)
    return zeros(nx, ny)
end

function jacobian(::Identity, nx, ny, args...)
    return I(min(nx, ny))
end

function jacobian(jac::Invariant, nx, ny, args...)
    return jac.value
end

function jacobian(jac::Constant, nx, ny, p)
    return jac.func(p)
end

function jacobian(jac::Constant, nx, ny, dx, x, p, t)
    return jac.func(p)
end

function jacobian(jac::Constant, nx, ny, dx, x, y, p, t)
    return jac.func(p)
end

function jacobian(jac::Union{Linear, Nonlinear}, nx, ny, args...)
    return jac.func(args...)
end

# in-place

function jacobian!(J, ::Empty, args...)
    return J
end

function jacobian!(J, ::Zeros, args...)
    return J .= 0
end

function jacobian!(J, ::Identity, args...)
    J .= 0
    for i = 1:size(J,1)
        J[i,i] = 1
    end
    return J
end

function jacobian!(J, jac::Invariant, args...)
    return J .= jac.value
end

function jacobian!(J, jac::Constant, p)
    return jac.func(J, p)
end

function jacobian!(J, jac::Constant, dx, x, p, t)
    return jac.func(J, p)
end

function jacobian!(J, jac::Constant, dx, x, y, p, t)
    return jac.func(J, p)
end

function jacobian!(J, jac::Union{Linear, Nonlinear}, args...)
    return jac.func(J, args...)
end

# --- Gradient Evaluation --- #

# out-of-place, statically sized

function gradient(::Empty, ::Val{NX}, args...) where NX
    return SVector{0, Float64}()
end

function gradient(::Zeros, ::Val{NX}, args...) where NX
    return zeros(SVector{NX, Float64})
end

function gradient(jac::Invariant, ::Val{NX}, args...) where NX
    return SVector{NX}(jac.value)
end

function gradient(jac::Constant, ::Val{NX}, p) where NX
    return SVector{NX}(jac.func(p))
end

function gradient(jac::Constant, ::Val{NX}, dx, x, p, t) where NX
    return SVector{NX}(jac.func(p))
end

function gradient(jac::Constant, ::Val{NX}, dx, x, y, p, t) where NX
    return SVector{NX}(jac.func(p))
end

function gradient(jac::Union{Linear, Nonlinear}, ::Val{NX}, args...) where NX
    return SVector{NX}(jac.func(args...))
end

# out-of-place, dynamically sized

function gradient(::Empty, nx, ny, args...)
    return Matrix{Float64}(val(nx), val(ny))
end

function gradient(::Zeros, nx, ny, args...)
    return zeros(val(nx), val(ny))
end

function gradient(::Identity, nx, ny, args...)
    return I(min(val(nx), val(ny)))
end

function gradient(grad::Invariant, nx, ny, args...)
    return grad.value
end

function gradient(grad::Constant, nx, ny, p)
    return grad.func(p)
end

function gradient(grad::Constant, nx, ny, dx, x, p, t)
    return grad.func(p)
end

function gradient(grad::Constant, nx, ny, dx, x, y, p, t)
    return grad.func(p)
end

function gradient(grad::Union{Linear, Nonlinear}, nx, ny, args...)
    return grad.func(args...)
end

# in-place

function gradient!(J, ::Empty, args...)
    return J
end

function gradient!(J, ::Zeros, args...)
    return J .= 0
end

function gradient!(J, ::Identity, args...)
    J .= 0
    for i = 1:size(J,1)
        J[i] = 1
    end
    return J
end

function gradient!(J, grad::Invariant, args...)
    return J .= grad.grad
end

function gradient!(J, grad::Constant, p)
    return grad.func(J, p)
end

function gradient!(J, grad::Constant, dx, x, p, t)
    return grad.func(J, p)
end

function gradient!(J, grad::Constant, dx, x, y, p, t)
    return grad.func(J, p)
end

function gradient!(J, grad::Union{Linear, Nonlinear}, args...)
    return grad.func(J, args...)
end
