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
