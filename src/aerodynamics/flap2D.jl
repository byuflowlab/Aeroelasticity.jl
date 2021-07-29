"""
    Flap2D <: NoStateModel

Linear, steady-state, two-dimensional control surface model with parameters
``c_{l,\delta}``, ``c_{d,\delta}``, and ``c_{m,\delta}.
"""
struct Flap2D <: NoStateModel end

"""
    Flap2D()

Initialize an object of type [`Flap2D`](@ref)
"""
Flap2D()

# --- Traits --- #

number_of_parameters(::Type{Flap2D}) = 2
