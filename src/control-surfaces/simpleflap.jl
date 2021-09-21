"""
    SimpleFlap <: NoStateModel

Linear, steady-state, two-dimensional control surface model with parameters
``c_{n,\\delta}``, ``c_{a,\\delta}``, and ``c_{m,\\delta}``.
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

function set_parameters!(p, model::SimpleFlap; cnd, cad, cmd)

    p[1] = cnd
    p[2] = cad
    p[3] = cmd

    return p
end

function separate_parameters(model::SimpleFlap, p)

    return (cnd=p[1], cad=p[2], cmd=p[3])
end

# --- Internal Methods --- #

function simpleflap_loads(b, U, ρ, cnδ, caδ, cmδ, δ)
    L = ρ*U^2*b*cnδ*δ
    D = ρ*U^2*b*caδ*δ
    M = 2*ρ*U^2*b^2*cmδ*δ
    return L, D, M
end
