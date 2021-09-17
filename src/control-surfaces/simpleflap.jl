"""
    SimpleFlap <: NoStateModel

Linear, steady-state, two-dimensional control surface model with parameters
``c_{l,\\delta}``, ``c_{d,\\delta}``, and ``c_{m,\\delta}``.
"""
struct SimpleFlap <: NoStateModel end

"""
    SimpleFlap()

Initialize an object of type [`SimpleFlap`](@ref)
"""
SimpleFlap()

# --- Traits --- #

number_of_parameters(::Type{SimpleFlap}) = 3

# --- Convenience Methods --- #

function set_parameters!(p, model::SimpleFlap; cld, cdd, cmd)

    p[1] = cld
    p[2] = cdd
    p[3] = cmd

    return p
end

function separate_parameters(model::SimpleFlap, p)

    return (cld=p[1], cdd=p[2], cmd=p[3])
end

# --- Internal Methods --- #

function simpleflap_loads(b, U, ρ, clδ, cdδ, cmδ, δ)
    L = ρ*U^2*b*clδ*δ
    D = ρ*U^2*b*cdδ*δ
    M = 2*ρ*U^2*b^2*cmδ*δ
    return L, D, M
end
