abstract type InPlaceness end
struct InPlace end
struct OutOfPlace end

function inplaceness(::Type{T}) where T
    if isinplace(T)
        return InPlace()
    else
        return OutOfPlace()
    end
end
inplaceness(model) = inplaceness(typeof(model))


abstract type MatrixState end
struct Identity end
struct Constant end
struct Varying end

abstract type InputDependence end
struct Linear end
struct Nonlinear end


function mass_matrix_state(::Type{T}) where T
    if has_mass_matrix(T)
        if constant_mass_matrix(T)
            return Constant()
        else
            return Varying()
        end
    else
        return Identity()
    end
end
mass_matrix_state(model) = mass_matrix_state(typeof(model))

function input_jacobian_state(::Type{T}) where T
    if constant_input_jacobian(T)
        return Constant()
    else
        return Varying()
    end
end
input_jacobian_state(model) = input_jacobian_state(typeof(model))

constant_matrix_state(::Identity) = true
constant_matrix_state(::Constant) = true
constant_matrix_state(::Varying) = false

function input_dependence(::Type{T}) where T
    if linear_input_dependence(T)
        return Linear()
    else
        return Nonlinear()
    end
end
input_dependence(model) = input_dependence(typeof(model))
