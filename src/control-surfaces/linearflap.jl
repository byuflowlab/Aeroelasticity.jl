"""
    LinearFlap <: NoStateModel

Linear, steady-state, two-dimensional control surface model with parameters
``c_{l,\\delta}``, ``c_{d,\\delta}``, and ``c_{m,\\delta}``.
"""
struct LinearFlap <: NoStateModel end

"""
    LinearFlap()

Initialize an object of type [`LinearFlap`](@ref)
"""
LinearFlap()

# --- Traits --- #

number_of_parameters(::Type{LinearFlap}) = 3

# --- Convenience Methods --- #

function set_parameters!(p, model::LinearFlap; cld, cdd, cmd)

    p[1] = cld
    p[2] = cdd
    p[3] = cmd

    return p
end

function separate_parameters(model::LinearFlap, p)

    return (cld=p[1], cdd=p[2], cmd=p[3])
end
