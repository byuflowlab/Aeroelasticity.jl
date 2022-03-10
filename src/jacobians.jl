# --- Jacobian Types --- #

abstract type JacobianType end

struct Empty <: JacobianType end
struct Zeros <: JacobianType end
struct Identity <: JacobianType end
struct Invariant{T} <: JacobianType; value::T; end
struct Constant{T} <: JacobianType; func::T; end
struct Linear{T} <: JacobianType; func::T; end
struct Nonlinear{T} <: JacobianType; func::T; end

# default constructors set jacobian (function) later
Invariant() = Invariant(nothing)
Constant() = Constant(nothing)
Linear() = Linear(nothing)
Nonlinear() = Nonlinear(nothing)

# trait dispatch
isempty(::JacobianType) = false
isempty(::Empty) = true

iszero(::JacobianType) = false
iszero(::Zeros) = true

isidentity(::JacobianType) = false
isidentity(::Identity) = true

isinvariant(::JacobianType) = false
isinvariant(::Empty) = true
isinvariant(::Zeros) = true
isinvariant(::Identity) = true
isinvariant(::Invariant) = true

isconstant(::JacobianType) = false
isconstant(::Empty) = true
isconstant(::Zeros) = true
isconstant(::Identity) = true
isconstant(::Invariant) = true
isconstant(::Constant) = true

islinear(::JacobianType) = false
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
