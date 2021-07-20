## Getting Started

```@setup guide
# this is placed here to pre-install matplotlib so the documentation doesn't get cluttered with the installation print statements.
using Plots
pyplot()
nothing #hide
```

In this guide we introduce you to the basic functionality of this package in a step by step manner.  This is a good starting point for learning about how to use this package.  For more details about how to use a particular model, it's [relevant documentation](@ref) is likely a better resource.  For more examples of how to use this package see the [examples](@ref Examples).

If you haven't yet, now would be a good time to install AerostructuralDynamics.  AerostructuralDynamics can be installed from the Julia REPL by typing `]` (to enter the package manager) and then running the following command.
```julia
pkg> add https://flow.byu.edu/AerostructuralDynamics.jl
```

Now, that the package is installed we need to load it so that we can use it.  It's also often helpful to load the LinearAlgebra package.

```@example guide
using AerostructuralDynamics, LinearAlgebra
nothing #hide
```

The geometry we will be working with is a typical section model with two degrees of freedom.
![](typical-section.svg)
The equations of motion for this model are
```math
\begin{bmatrix} m & S_\theta \\ S_\theta & I_\theta \end{bmatrix}
\begin{Bmatrix} \ddot{h} \\ \ddot{\theta} \end{Bmatrix} +
\begin{bmatrix} K_h & 0 \\ 0 & K_\theta \end{bmatrix}
\begin{Bmatrix} h \\ \theta \end{Bmatrix} =
\begin{Bmatrix} -\mathcal{L} \\ \mathcal{M} \end{BMatrix}
```
where ``a`` is the normalized location of the reference point relative to the semi-chord, ``b`` is the semi-chord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the static imbalance, ``I_θ`` is the moment of inertia about the reference point, ``\mathcal{L}`` is the lift per unit span, and ``\mathcal{M}`` is the moment per unit span about the reference point.

## Initializing a Model

To simulate the behavior of the typical section model, we will be using the [`TypicalSection`](@ref) model.  We can obtain details about this model from [its documentation](@ref Typical Section Model).

```@example guide
model = TypicalSection()
```

This model is an explicit model, which means that its governing equations may be expressed as
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x,y,p,t)`` is a function which defines the mass matrix corresponding to the differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of states, ``y`` is a vector of inputs/coupling variables, ``p`` is a vector of parameters, and ``t`` is the current time.  If this model was an implicit model, its governing equations would be expressed as
```math
0 = f(\dot{x},x,y,p,t)
```

# Defining State Variables, Inputs, and/or Parameters

As noted [in its documentation](@ref Typical Section Model), the state variable, input, and parameter vectors for the [`TypicalSection`](@ref) model are defined as
```math
u = \begin{bmatrix} h \\ \theta \\ \dot{h} \\ \dot{\theta} \end{bmatrix} \quad y = \begin{bmatrix} \mathcal{L} \\ \mathcal{M} \end{bmatrix} \quad p = \begin{bmatrix} k_h \\ k_\theta \\ m \\ S_\theta \\ I_\theta \end{bmatrix}
```

For this guide, we define these vectors as follows:

```@example guide
# initial state variable vector
h = 0
θ = 0
hdot = 0
θdot = 0
x = [h, θ, hdot, θdot]

# initial input vector
L = 0
M = 0
y = [L, M]

# initial parameter vector
kh = 1
kθ = 1
m = 1
Sθ = 1
Iθ = 1
p = [kh, kθ, m, Sθ, Iθ]

nothing #hide
```

In this case the small number of states, inputs, and parameters make manually constructing the state, input, and parameter vectors easy.  However, for models with larger numbers of states, inputs, and/or parameters it is often more convenient to define the values contained in these vectors in a different format.  In this case, specialized versions of the
[`to_state_vector`](@ref),  [`to_input_vector`](@ref), [`to_parameter_vector`](@ref) [`from_state_vector`](@ref),  [`from_inputs`](@ref), and/or [`from_parameters`](@ref) functions for each model may be used (or their equivalent in-place versions) to convert to/from the vectorial representation of each model's states, inputs, and/or parameters.

```@example guide
x = to_state_vector(model, h, θ, hdot, θdot)
y = to_input_vector(model, L, M)
p = to_parameter_vector(model, kh, kθ, m, Sθ, Iθ)
```

# Performing a Stability Analysis

The stability of a model for a given set of inputs and parameters may be determined by calling the [`stability_analysis`](@ref) function.  For nonlinear systems, the provided state variables must correspond to an equilibrium point for the stability analysis to be theoretically valid.

```@example guide
λ, V = stability_analysis(model, u, y, p)
```

# Performing a Simulation

A model may be simulated by first defining an [`DifferentialEquations.ODEProblem`](@ref) (or [`DifferentialEquations.DAEProblem`](@ref) if the model is implicit) which may then be solved using the [`DifferentialEquations`](@ref) package.  To facilitate this process, this package provides convenience constructors for the [`DifferentialEquations.ODEProblem`](@ref) and [`DifferentialEquations.DAEProblem`](@ref) functions.

```@example guide
using DifferentialEquations

prob = DifferentialEquations.ODEProblem(model, u0, tspan, y, p; kwargs...)

sol = DifferentialEquations.solve(prob)
```
