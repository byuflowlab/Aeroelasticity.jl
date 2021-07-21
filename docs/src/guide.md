# Getting Started

```@setup guide
# this is placed here to pre-install matplotlib so the documentation doesn't get cluttered with the installation print statements.
using Plots
pyplot()
nothing #hide
```

In this guide we introduce you to the basic functionality of this package in a step by step manner.  This is a good starting point for learning about how to use this package.  For more details about how to use a particular model, the model documentation is likely a better resource.  For more examples of how to use this package see the [examples](@ref Examples).

If you haven't yet, now would be a good time to install AerostructuralDynamics.  AerostructuralDynamics can be installed from the Julia REPL by typing `]` (to enter the package manager) and then running the following command.
```julia
pkg> add https://flow.byu.edu/AerostructuralDynamics.jl
```

Now, that the package is installed we need to load it so that we can use it.

```@example guide
using AerostructuralDynamics
nothing #hide
```

The geometry we will be working with is a typical section model with two degrees of freedom.

![](typical-section.svg)

We will be coupling this model with an unsteady aerodynamic model in order to create and analyze a 2D aeroelastic system.

## Initializing Models

We need to initialize two models to create our aeroelastic system: Peter's finite state aerodynamic model [`Peters`](@ref) and a typical section structural model [`TypicalSection`](@ref).  A detailed description of these models may be found in the [model documentation](@ref Model Documentation).  We couple these models together by calling the [`couple_models`](@ref) function.  

```@example guide
aerodynamic_model = Peters{4}() # we use four aerodynamic state variables
structural_model = TypicalSection()
coupled_model = couple_models(aerodynamic_model, structural_model)
nothing #hide
```

The governing equations for both of these models may be expressed as the first order ordinary differential equation
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x,y,p,t)`` is a function which defines the mass matrix corresponding to the differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of states, ``y`` is a vector of inputs/coupling variables, ``p`` is a vector of parameters, and ``t`` is the current time.

## Defining State Variables, Inputs, and Parameters

As described in [its documentation](@ref peters-finite-state-model), the state, input, and parameter vectors for the aerodynamic model are
```math
x_\text{aero} = \begin{bmatrix} \lambda_1 \\ \lambda_2 \\ \vdots \\ \lambda_N \end{bmatrix} \quad
y_\text{aero} = \begin{bmatrix} u \\ v \\ \omega \end{bmatrix} \quad
p = \begin{bmatrix} a \\ b \\ a_0 \\ \alpha_0 \end{bmatrix}
```
where ``\lambda_1, \lambda_2, \dots, \lambda_N`` are the aerodynamic states,
``u`` is the chordwise freestream velocity, ``v`` is the normal freestream velocity, ``\omega`` is the angular freestream velocity, ``a`` is the normalized reference location relative to the semi-chord, ``b`` is the semi-chord, ``a_0`` is the section lift slope, and ``\alpha_0`` is the section zero lift angle of attack.  Positive freestream velocity components are defined as shown in the following figure.
![](airfoil.svg)

As described in [its documentation](@ref typical-section-model), the state, input, and parameter vectors for the typical section model are defined as
```math
x_\text{stru} = \begin{bmatrix} h \\ \theta \\ \dot{h} \\ \dot{\theta} \end{bmatrix} \quad y_\text{stru} = \begin{bmatrix} \mathcal{L} \\ \mathcal{M} \end{bmatrix} \quad p_\text{stru} = \begin{bmatrix} k_h \\ k_\theta \\ m \\ S_\theta \\ I_\theta \end{bmatrix}
```
where ``h`` is plunge, ``\theta`` is pitch, ``\mathcal{L}`` is the lift per unit span, ``\mathcal{M}`` is the moment per unit span about the reference point, ``k_h`` is the linear spring constant, and ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the structural imbalance, and ``I_θ`` is the mass moment of inertia about the reference point.

The state variables and inputs of a coupled model correspond to the state variables and inputs of its component models concatenated.  The parameters of a coupled model correspond to the parameters of its component models concatenated, followed by additional parameters which are specific to the coupled models.  The additional parameters introduced by the coupled model in this case is the freestream velocity ``U`` and air density ``\\rho``.  The state, input, and parameter vectors for the coupled model are therefore
```math
x_\text{coupled} = \begin{bmatrix} \lambda_1 \\ \lambda_2 \\ \vdots \\ \lambda_N \\ h \\ \theta \\ \dot{h} \\ \dot{\theta} \end{bmatrix} \quad y_\text{coupled} = \begin{bmatrix} u \\ v \\ \omega \\ \mathcal{L} \\ \mathcal{M} \end{bmatrix} \quad p_\text{coupled} = \begin{bmatrix} a \\ b \\ a_0 \\ \alpha_0 \\ k_h \\ k_\theta \\ m \\ S_\theta \\ I_\theta \\ U \\ \rho \end{bmatrix}
```

For this guide, we use the following state variables and parameters
```@example guide
# state variables
x_coupled = zeros(number_of_states(coupled_model))

# non-dimensional parameters
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
V = 1.0 # = U/(b*ωθ) (reduced velocity)

# chosen dimensional parameters
b = 1
ρ = 1
ωθ = 1

# derived dimensional parameters
m = μ*ρ*pi*b^2
Sθ = m*xθ*b
Iθ = r2*m*b^2
ωh = σ*ωθ
kh = m*ωh^2
kθ = Iθ*ωθ^2
U = V*b*ωθ

# parameter vector
p_aero = [a, b, a0, α0]
p_stru = [kh, kθ, m, Sθ, Iθ]
p_additional = [U, ρ]
p_coupled = vcat(p_aero, p_stru, p_additional)

nothing #hide
```

For coupled models, inputs are calculated as a function of the coupled model's states and parameters, as well as the current time.  For uncoupled models, inputs must be manually specified.

```@example guide
t = 0 # time

y_coupled = get_inputs(coupled_model, x_coupled, p_coupled, t) # inputs
nothing #hide
```

For this model, the small number of states, inputs, and parameters make manually constructing and inspecting the state, input, and parameter vectors easy.  However, for models with larger numbers of states, inputs, and/or parameters, converting to/from a vectorial representation of the states, inputs, and/or parameters may be challenging.  For this reason, specialized versions of the [`to_state_vector`](@ref),  [`to_input_vector`](@ref), [`to_parameter_vector`](@ref) [`from_state_vector`](@ref),  [`from_inputs`](@ref), and/or [`from_parameters`](@ref) functions for each model and/or model coupling may be used (or their equivalent in-place versions) to convert to/from the vectorial representation of each model's states, inputs, and/or parameters.

## Performing a Stability Analysis

The stability of a model for a given set of state variables, inputs, and parameters may be determined by calling the [`stability_analysis`](@ref) function, which determines .  For nonlinear systems, the provided state variables must correspond to an equilibrium point for the stability analysis to be theoretically valid.  In our case, our aeroelastic system is linear, so no such need exists.

```@example guide
λ, U, V = stability_analysis(coupled_model, x_coupled, y_coupled, p_coupled, t)
nothing #hide
```

A positive real part corresponding to any eigenvalue returned from the [`stability_analysis`](@ref) function would indicate that the system is unstable.

## Performing a Simulation

To simulate the behavior of our model we first need to create an object of type [`DifferentialEquations.ODEFunction`](@ref) using [`ode_function`](@ref).  Then the [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl) package may be used to solve the ordinary differential equation corresponding to the model.

```@example guide
using DifferentialEquations

# non-zero plunge degree of freedom
u0 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

# simulate for 10 seconds
tspan = (0.0, 10.0)

# construct ODE function
f = ode_function(coupled_model)

# construct ODE problem
prob = DifferentialEquations.ODEProblem(f, u0, tspan, p_coupled)

# solve ODE
sol = DifferentialEquations.solve(prob)

nothing #hide
```

## Visualizing Results
