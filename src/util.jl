# --- extract values from Val() objects --- #

val(x) = x
val(::Val{x}) where x = x

# --- automatic differentiation --- #

# out-of-place
function autodiff_jacobian_func(f, ivar)
    (args...) -> begin
        x = args[ivar]
        fy = (x) -> f(Base.setindex(args, x, ivar)...)
        ForwardDiff.jacobian(fy, x)
    end
end

# in-place
function autodiff_jacobian_func(f, ny, ivar)
    (J, args...) -> begin
        x = args[ivar]
        y = similar(J, val(ny))
        fy = (y, x) -> f(y, Base.setindex(args, x, ivar)...)
        ForwardDiff.jacobian!(J, fy, y, x)
    end
end

# out-of-place
function autodiff_derivative_func(f, ivar)
    (args...) -> begin
        x = args[ivar]
        fy = (x) -> f(Base.setindex(args, x, ivar)...)
        ForwardDiff.derivative(fy, x)
    end
end

# in-place
function autodiff_derivative_func(f, ny, ivar)
    (dy, args...) -> begin
        x = args[ivar]
        y = similar(dy, val(ny))
        fy = (y, x) -> f(y, Base.setindex(args, x, ivar)...)
        ForwardDiff.derivative!(dy, fy, y, x)
    end
end

# --- block diagonal --- #

function block_diagonal(Jii)
    
    J = zeros(sum(size.(Jii, 1)), sum(size.(Jii, 2)))

    i = j = 1
    for iblock = 1:N
        nx, ny = size(Jii[iblock])
        ix, iy = i:i+nx-1, j:j+ny-1
        J[ix, iy] = Jii[iblock]
        i += nx; j += ny
    end
    
    return J
end

block_diagonal(Jii::NTuple{N,StaticMatrix}) where N = static_block_diagonal(Jii)

@generated function static_block_diagonal(Jii::NTuple{N,StaticMatrix}) where N

    # assemble first row of the jacobian matrix
    expr = quote
        J = Jii[1]
    end
    for j = 2:N
        if 1 == j
            expr = quote
                $expr
                Jij = Jii[1]
            end
        else
            expr = quote
                $expr
                nx = Size(Jii[1])[1]
                ny = Size(Jii[$j])[2]
                sz = Size(nx, ny)
                Jij = zero(similar_type(Jii[1], sz))
            end
        end
        expr = quote
            $expr
            J = hcat(J, Jij)
        end
    end

    # add remaining rows of the jacobian matrix
    for i = 2:N
        if i == 1
            expr = quote
                $expr
                Ji = Jii[$i]
            end
        else
            expr = quote
                $expr
                nx = Size(Jii[$i])[1]
                ny = Size(Jii[1])[2]
                sz = Size(nx, ny)
                Ji = zero(similar_type(Jii[$i], sz))
            end
        end
        for j = 2:N
            if i == j
                expr = quote
                    $expr
                    Jij = Jii[$i]
                end
            else
                expr = quote
                    $expr
                    nx = Size(Jii[$i])[1]
                    ny = Size(Jii[$j])[2]
                    sz = Size(nx, ny)
                    Jij = zero(similar_type(Jii[$i], sz))
                end
            end
            expr = quote
                $expr
                Ji = hcat(Ji, Jij)
            end
        end
        expr = quote
            $expr
            J = vcat(J, Ji)
        end
    end

    return expr
end
