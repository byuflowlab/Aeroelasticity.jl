# Getting Started

```@setup guide
# this is placed here to pre-install matplotlib so the documentation doesn't get cluttered with the installation print statements.
using Plots
pyplot()
nothing #hide
```

In this guide we introduce you to the basic functionality of this package in a step by step manner.  This is a good starting point for learning about how to use this package.  For more details about how to use a particular model or combination models, refer to the model documentation.  For more examples of how to use this package see the [examples](@ref Examples).

```@contents
Pages = ["library.md"]
Depth = 3
```

## Overview

AerostructuralDynamics is designed to simulate and assess the stability of complex aerostructural systems.  It does this by constructing a system of first order ordinary differential equations and associated jacobians for various aerostructural systems which may be used to find equilibrium points, perform stability analyses, and/or perform time domain simulations.

This package is designed to be modular, so that models may be easily swapped out for alternative models.  It is also designed to be extensible, so that new models may be incorporated with minimal effort.  To achieve these aims, this package assumes that the governing equations for all models may be expressed as the implicit ordinary differential equation
```math
0 = f(\dot{x},x,y,p,t)
```
where ``f(\dot{x}, x, y, p, t)`` is a residual function, ``\dot{x}`` is a vector of state rates, ``x`` is a vector of state variables, ``y`` is a vector of time-varying parameters, ``p`` is a vector of time-independent parameters, and ``t`` is the current time.  For the purposes of this package we will refer to the time-varying parameters ``y`` as "inputs" and the time-independent parameters ``p`` as "parameters".

Since the governing equations for all models in this package follow the same general form, the process for performing analyses using any standalone and/or coupled model provided by this package is the same.  First, the relevant standalone and/or coupled model must be initialized.  Then, the model's initial state rates, states, inputs, parameters, and time must be defined.  Finally, the chosen analysis may be performed.  We demonstrate how to perform these steps in the following sections.

## Initializing a Model

If you haven't yet, now would be a good time to install AerostructuralDynamics.  It can be installed from the Julia REPL by typing `]` (to enter the package manager) and then running the following command.
```julia
pkg> add https://flow.byu.edu/AerostructuralDynamics.jl
```

Now, that the package is installed we need to load it so that we can use it.

```@example guide
using AerostructuralDynamics
nothing #hide
```

For the purposes of this guide, we will be working with a two-degree-of-freedom typical section model, as shown in the following figure.

![](typical-section.svg)

Our goal is to create an 2D aeroelastic model which we can use to simulate the behavior of this system.

For the aerodynamic model, we will be using Peters' finite state model with four aerodynamic state variables (see [`Peters`](@ref)).  For the structural model, we will be using the typical section model (see [`TypicalSection`](@ref)).  To create a coupled model using these model, we use the [`couple_models`](@ref) function.  Details about how to initialize these models may be found in the model documentation.

```@example guide
# initialize the aerodynamic model
aerodynamic_model = Peters{4}()

# initialize the structural model
structural_model = TypicalSection()

# initialize the structural model
coupled_model = couple_models(aerodynamic_model, structural_model)

nothing #hide
```

## Defining State Rates, State Variables, Inputs, and Parameters

As described in the documentation for the [`Peters`](@ref) model, its state, input, and parameter vectors are defined as
```math
x_\text{aero} = \begin{bmatrix} \lambda_1 \\ \lambda_2 \\ \vdots \\ \lambda_N \end{bmatrix} \quad
y_\text{aero} = \begin{bmatrix} u \\ v \\ \omega \end{bmatrix} \quad
p = \begin{bmatrix} a \\ b \\ a_0 \\ \alpha_0 \end{bmatrix}
```
where ``\lambda_1, \lambda_2, \dots, \lambda_N`` are the aerodynamic states,
``u`` is the chordwise freestream velocity, ``v`` is the normal freestream velocity, ``\omega`` is the angular freestream velocity, ``a`` is the normalized reference location relative to the semi-chord, ``b`` is the semi-chord, ``a_0`` is the section lift slope, and ``\alpha_0`` is the section zero lift angle of attack.  Positive freestream velocity components are defined as shown in the following figure.

![](airfoil.svg)

As described in the documentation for the [`TypicalSection`](@ref) model, the state, input, and parameter vectors for the typical section model are defined as
```math
x_\text{stru} = \begin{bmatrix} h \\ \theta \\ \dot{h} \\ \dot{\theta} \end{bmatrix} \quad y_\text{stru} = \begin{bmatrix} \mathcal{L} \\ \mathcal{M} \end{bmatrix} \quad p_\text{stru} = \begin{bmatrix} k_h \\ k_\theta \\ m \\ S_\theta \\ I_\theta \end{bmatrix}
```
where ``h`` is plunge, ``\theta`` is pitch, ``\mathcal{L}`` is the lift per unit span, ``\mathcal{M}`` is the moment per unit span about the reference point, ``k_h`` is the linear spring constant, and ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the structural imbalance, and ``I_θ`` is the mass moment of inertia about the reference point.

In this package, we define the state variables and inputs of a coupled model as the state variables and inputs of its component models concatenated.  We also define the parameters of a coupled model as the parameters of its component models concatenated, followed by a set of additional parameters which are specific to the coupled model.  As noted in the documentation for the coupled model we are using in this guide, the additional parameters introduced by the coupled model we consider in this example are the freestream velocity ``U_\infty`` and air density ``\rho``.  The state, input, and parameter vectors for this coupled model are therefore
```math
x_\text{coupled} = \begin{bmatrix} \lambda_1 \\ \lambda_2 \\ \vdots \\ \lambda_N \\ h \\ \theta \\ \dot{h} \\ \dot{\theta} \end{bmatrix} \quad y_\text{coupled} = \begin{bmatrix} u \\ v \\ \omega \\ \mathcal{L} \\ \mathcal{M} \end{bmatrix} \quad p_\text{coupled} = \begin{bmatrix} a \\ b \\ a_0 \\ \alpha_0 \\ k_h \\ k_\theta \\ m \\ S_\theta \\ I_\theta \\ U_\infty \\ \rho \end{bmatrix}
```

For standalone models, both inputs ``y`` parameters ``p`` are user-specified.  For coupled models, however, inputs to one model in the coupled system often correspond to outputs from another model in the coupled system.  For example, the lift and moment required by the typical section model are defined by Peter's finite state model.    correspond to 

  For coupled models, inputs ``y`` are defined as a function of the state rates, states, and parameters of all the models in the coupled system and the current time.  

 function of the state rates, states, and Parameters
only parameters ``p`` are user-specified while  inputs are defined as a function of the  assumes that the inputs to each of its subcomponent models may be defined as a function of the  

however, inputs to a model can often be defined as a function of the state variables and parameters of the other models in the coupled system.  The interdependencies between models may therefore be captured by allowing the inputs to each model to be defined as a function of the state variables and parameters of the coupled system.  Each model in a coupled model is therefore coupled together through their inputs, which

Rate, state, input, and parameter vectors may either be constructed directly or initialized using the [`get_states`](@ref), [`get_inputs`](@ref), and/or [`get_parameters`](@ref) convenience functions.  For coupled models, only the state and parameter vectors must be defined since the inputs for coupled models are defined as a function of the coupled model's state variables, parameters, and current time.

```@example guide
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

# initial rates
dx_coupled = zeros(number_of_states(coupled_model))

# initial states
x_coupled = zeros(number_of_states(coupled_model))

# initial parameters
p_aero = [a, b, a0, α0]
p_stru = [kh, kθ, m, Sθ, Iθ]
p_additional = [U, ρ]
p_coupled = vcat(p_aero, p_stru, p_additional)

# initial time
t = 0

nothing #hide
```

To calculate the inputs corresponding to a given set of state variables and parameters for a coupled model the [`get_coupling_inputs`](@ref) function may be used.

```@example guide
# calculate input vector
y_coupled = get_coupling_inputs(coupled_model, dx_coupled, x_coupled, p_coupled, t) # inputs

nothing #hide
```

In order to more easily interpret the elements of the rate, state, input, and/or parameter vectors, the [`separate_states`](@ref), [`separate_inputs`](@ref), and/or [`separate_parameters`](@ref) functions may be used.  These functions separate and assign names to the elements of the state, input and/or parameter vectors so that the identity of each element in these vectors may be more easily understood.

```@example guide
rates = separate_states(coupled_model, dx_coupled)
states = separate_states(coupled_model, x_coupled)
inputs = separate_inputs(coupled_model, y_coupled)
parameters = separate_parameters(coupled_model, p_coupled)

nothing #hide
```





  defined as a function of the state variables and parameters of the coupled system.

for each model often correspond to outputs from other models of the state variables and parameters of the coupled system.


When models are coupled together, this package defines the inputs ``y`` using output functions corresponding to other models.  For example, applied loads may be considered an input to a structural

to each model as a function of the outputs from other models.
In couple models together, this package assumes that the inputs ``y`` to each model correspond to outputs from other models.  The outputs  which may be defined as a function of the state rates, states, and parameters of any and/or all of the coupled models as well as the current time.





## Performing a Stability Analysis

The stability of a model for a given set of state variables, inputs, and parameters may be determined by calling the [`get_eigen`](@ref) function, which returns eigenvalues, left eigenvectors, and right eigenvectors.  For nonlinear systems, the provided state variables must correspond to an equilibrium point for the stability analysis to be theoretically valid.  In our case, our aeroelastic system is linear with respect to the state variables, so any set of state variables will yield the same result.

```@example guide
λ, U, V = get_eigen(coupled_model, dx_coupled, x_coupled, y_coupled, p_coupled, t)
nothing #hide
```

A positive real part corresponding to any eigenvalue returned from the [`get_eigen`](@ref) function indicates that the system is unstable for the provided set of state variables, inputs, and parameters.

## Performing a Simulation

To simulate the behavior of our model we first need to create an object of type `DifferentialEquations.ODEFunction` using the [`get_ode`](@ref) function.  Then the [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl) package may be used to solve the ordinary differential equation corresponding to the model.

```@example guide
using DifferentialEquations

# non-zero plunge degree of freedom
u0 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

# simulate for 10 seconds
tspan = (0.0, 100.0)

# construct ODE function
f = get_ode(coupled_model)

# construct ODE problem
prob = DifferentialEquations.ODEProblem(f, u0, tspan, p_coupled)

# solve ODE
sol = DifferentialEquations.solve(prob)

```

# Visualizing Simulation Results

We can use DifferentialEquations' built-in interface with the [Plots](https://github.com/JuliaPlots/Plots.jl) package to plot the simulation results.

```@example guide
using Plots
pyplot()

plot(sol,
    vars = [5,6,7,8],
    xlabel = "t",
    ylabel = permutedims([
        "\$h\$",
        "\$\\theta\$",
        "\$\\dot{h}\$",
        "\$\\dot{\\theta}\$",
        ]),
    label = "",
    layout = (4, 1),
    size = (600,1200)
    )

savefig("guide-solution.svg") #hide

nothing #hide
```

![](guide-solution.svg)

For some models, we can also visualize the solution geometry and/or create animations with the help of custom plot recipes provided by this package.

```@example typical-section-stability
# Plot Recipes:
# - plot(model, dx, x, p, t)
# - plot(model, sol) # plots the solution geometry at index `sol.tslocation`
# - plot(model, sol, t) # plots the solution geometry at time `t`

# create animation
anim = @animate for t in range(tspan[1], tspan[2], length=200)
    plot(coupled_model, sol, t)
end

# save animation
gif(anim, "guide-simulation.gif")

nothing #hide
```

![](guide-simulation.gif)
