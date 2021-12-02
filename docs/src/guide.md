# Getting Started

```@setup guide
# this is placed here to pre-install matplotlib so the documentation doesn't get cluttered with the installation print statements.
using Plots
pyplot()
nothing #hide
```

In this guide we introduce you to the basic functionality of this package in a step by step manner.  This is a good starting point for learning about how to use this package.  For more details about how to use a particular model and/or coupling, refer to the documentation for the model and/or coupling.  For more examples of how to use this package see the [examples](@ref Examples).

```@contents
Pages = ["library.md"]
Depth = 3
```

## Overview

AerostructuralDynamics is designed to simulate and assess the stability of complex aerostructural systems.  It does this by constructing a system of first order ordinary differential equations and associated jacobians for various aerostructural systems which may be used to find equilibrium points, perform stability analyses, and/or perform time domain simulations.

This package is designed to be modular, so that models may be easily swapped out for other models and/or customized.  It is also designed to be extensible, so that new models may be incorporated with minimal effort.  To achieve these aims, this package assumes that the governing equations for all models may be expressed as the implicit ordinary differential equation
```math
0 = f(\dot{x},x,y,p,t)
```
where ``f(\dot{x}, x, y, p, t)`` is a residual function, ``\dot{x}`` is a vector of state rates, ``x`` is a vector of state variables, ``y`` is a vector of time-varying parameters, ``p`` is a vector of time-independent parameters, and ``t`` is the current time.  For the purposes of this package we will refer to the time-varying parameters ``y`` as "inputs" and the time-independent parameters ``p`` as "parameters".

Since the governing equations for all models in this package follow the same general form, the process for performing analyses using any model provided by this package is the same.  First, the relevant model must be initialized.  Then, the model's initial state rates, states, inputs, parameters, and time are defined.  Finally, the chosen analysis is performed.  We demonstrate how to perform these steps in the following sections.

## Initializing Models

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

For the aerodynamic model, we will be using Peters' finite state model with four aerodynamic state variables (see [`Peters`](@ref)).  For the structural model, we will be using the typical section model (see [`TypicalSection`](@ref)).  The manner in which these two models will be coupled is defined by the [`PetersCoupling`](@ref) function.  To create a coupled model, we use the submodels and coupling to construct an object of type [`CoupledModel`](@ref).

```@example guide
# define the aerodynamic model
aerodynamic_model = Peters(N)

# define the structural model
structural_model = TypicalSection()

# define the submodels
submodels = (aerodynamic_model, structural_model)

# define the coupling between the two models
coupling = PetersCoupling(N)

# construct the coupled model
model = CoupledModel(submodels, coupling)

nothing #hide
```

## Defining State Rates, State Variables, Inputs, and Parameters

As described in the documentation for the [`Peters`](@ref) model, its state, input, and parameter vectors are defined as
```math
x_\text{aero} = \begin{bmatrix} \lambda_1 \\ \lambda_2 \\ \vdots \\ \lambda_N \end{bmatrix} \quad
y_\text{aero} = \begin{bmatrix} u \\ v \\ \omega \end{bmatrix} \quad
p = \begin{bmatrix} a \\ b \\ a_0 \\ \alpha_0 \\ c_{d_0} \\ c_{m_0} \end{bmatrix}
```
where ``\lambda_1, \lambda_2, \dots, \lambda_N`` are the aerodynamic states,
``u`` is the chordwise freestream velocity, ``v`` is the normal freestream velocity, ``\omega`` is the angular freestream velocity, ``a`` is the normalized reference location relative to the semi-chord, ``b`` is the semi-chord, ``a_0`` is the section lift slope, ``\alpha_0`` is the section zero lift angle of attack, ``c_{d_0}`` is the zero lift drag, and ``c_{m_0}`` is the zero lift moment.  Positive freestream velocity components are defined as shown in the following figure.

![](airfoil.svg)

As described in the documentation for the [`TypicalSection`](@ref) model, the state, input, and parameter vectors for the typical section model are defined as
```math
x_\text{stru} = \begin{bmatrix} h \\ \theta \\ \dot{h} \\ \dot{\theta} \end{bmatrix} \quad y_\text{stru} = \begin{bmatrix} \mathcal{L} \\ \mathcal{M} \end{bmatrix} \quad p_\text{stru} = \begin{bmatrix} k_h \\ k_\theta \\ m \\ S_\theta \\ I_\theta \end{bmatrix}
```
where ``h`` is plunge, ``\theta`` is pitch, ``\mathcal{L}`` is the lift per unit span, ``\mathcal{M}`` is the moment per unit span about the reference point, ``k_h`` is the linear spring constant, and ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the structural imbalance, and ``I_θ`` is the mass moment of inertia about the reference point.

In this package, we define the state variables and inputs of a coupled model as the state variables and inputs of its submodels concatenated.  We also define the parameters of a coupled model as the parameters of its submodels concatenated, followed by a set of additional parameters which are specific to the chosen coupling.  As noted in the documentation for [`PetersCoupling`](@ref), the additional parameters introduced by the coupled model we consider in this example are the freestream velocity ``U_\infty``, air density ``\rho``, and air speed of sound ``c``.  The state, input, and parameter vectors for this coupled model are therefore
```math
x_\text{coupled} = \begin{bmatrix} \lambda_1 \\ \lambda_2 \\ \vdots \\ \lambda_N \\ h \\ \theta \\ \dot{h} \\ \dot{\theta} \end{bmatrix} \quad y_\text{coupled} = \begin{bmatrix} u \\ v \\ \omega \\ \mathcal{L} \\ \mathcal{M} \end{bmatrix} \quad p_\text{coupled} = \begin{bmatrix} a \\ b \\ a_0 \\ \alpha_0 \\ c_{d_0} \\ c_{m_0} \\ k_h \\ k_\theta \\ m \\ S_\theta \\ I_\theta \\ U_\infty \\ \rho \\ c \end{bmatrix}
```

Rate, state, input, and parameter vectors may either be constructed directly or initialized using the [`get_states`](@ref), [`get_inputs`](@ref), and/or [`get_parameters`](@ref) convenience functions.  The in-place equivalents to these functions: [`set_states!`](@ref), [`set_inputs!`](@ref), and/or [`set_parameters!`](@ref) may also be used.  In the following block of code, we set the parameters for our coupled model directly.

```@example guide
# non-dimensional parameters
V = 1.0 # = U/(b*ωθ) (reduced velocity)
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a # normalized distance from the reference point to the center of mass
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
cd0 = 0 # zero lift drag coefficient
cm0 = 0 # zero lift moment coefficient

# chosen dimensional parameters
b = 1 # semi-chord
ρ = 1 # air density
ωθ = 1 # pitch natural frequency

# derived dimensional parameters
U = V*b*ωθ # velocity
m = μ*ρ*pi*b^2 # mass
Sθ = m*xθ*b # structural imbalance
Iθ = r2*m*b^2 # moment of inertia
ωh = σ*ωθ # plunge natural frequency
kh = m*ωh^2 # linear spring constant
kθ = Iθ*ωθ^2 # torsional spring constant

# parameters
p_aero = [a, b, a0, α0, cd0, cm0]
p_stru = [kh, kθ, m, Sθ, Iθ]
p_additional = [U, ρ]
p = vcat(p_aero, p_stru, p_additional)

nothing #hide
```

For standalone models, both inputs ``y`` parameters ``p`` are user-specified.  For coupled models, however, inputs to one model in the coupled system often correspond to outputs from another model in the coupled system.  For example, the lift and moment required by the typical section model are defined by Peter's finite state model and the velocities required by Peters' finite state model are defined as a function of the state variables of the typical section model.  To model these interdependencies, this package assumes that the inputs of a coupled model may be defined as a function of the state rates, states, and parameters of the coupled model as well as the current time.
```math
y = g(\dot{x}, x, p, t)
```
To evaluate this function to find the inputs for a coupled model, the [`get_coupling_inputs`](@ref) function may be used.

```@example guide
# choose rates, state, and time at which to calculate the coupling inputs
dx = zeros(number_of_states(model))
x = zeros(number_of_states(model))
t = 0

# calculate the coupling inputs
y = get_coupling_inputs(model, dx, x, p, t)

nothing #hide
```

In order to more easily interpret the elements of the rate, state, input, and/or parameter vectors, the [`separate_states`](@ref), [`separate_inputs`](@ref), and/or [`separate_parameters`](@ref) functions may be used.  These functions separate and assign names to the elements of the state, input and/or parameter vectors so that the identity of each element in these vectors may be more easily understood.

```@example guide
rates = separate_states(model, dx)
states = separate_states(model, x)
inputs = separate_inputs(model, y)
parameters = separate_parameters(model, p)

nothing #hide
```

## Finding an Equilibrium Point

To find equilibrium points, we first need to create an object of type `DifferentialEquations.ODEFunction` using the [`get_ode`](@ref) function.

```@example guide
# returns an ODEFunction
f = get_ode(model)

nothing #hide
```

Then a steady state solution may be found using DifferentialEquations.

```@example guide
using DifferentialEquations

# initial guess for state variables
x0 = zeros(number_of_states(model))

# steady state problem
prob = SteadyStateProblem(f, x0, p)

# steady state solution
x_ss = solve(prob, SSRootfind())

nothing #hide
```

Since our system is linear with respect to the state variables it has one equilibrium point at the origin.

Note that the coupling inputs are calculated automatically as a function of the state variables and parameters for coupled models.  For standalone model evaluation, the inputs are appended to the end of the parameter vector.

## Performing a Stability Analysis

The stability of a model for a given set of state variables, inputs, and parameters may be determined by calling the [`get_eigen`](@ref) function, which returns eigenvalues, left eigenvectors, and right eigenvectors.  For nonlinear systems, the provided state variables must correspond to an equilibrium point for the stability analysis to be theoretically valid.  Since our aeroelastic system is linear with respect to the state variables, any set of state variables will yield the same result.

```@example guide
λ, U, V = get_eigen(model, x_ss, p)

nothing #hide
```

A positive real part corresponding to any eigenvalue returned from the [`get_eigen`](@ref) function indicates that the system is unstable for the provided set of state variables, inputs, and parameters.


## Performing a Time Domain Simulation

To simulate the behavior of our model we first need to create an object of type `DifferentialEquations.ODEFunction` using the [`get_ode`](@ref) function.  

```@example guide
f = get_ode(model)

nothing #hide
```

Then the [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl) package may be used to solve the ordinary differential equation corresponding to the model.

```@example guide
using DifferentialEquations

# non-zero plunge degree of freedom
x0 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

# simulate for 100 seconds
tspan = (0.0, 100.0)

# construct ODE problem
prob = DifferentialEquations.ODEProblem(f, x0, tspan, p)

# solve ODE
sol = DifferentialEquations.solve(prob)

nothing #hide
```

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