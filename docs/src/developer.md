# Developer's Guide

In this guide, we describe how to create new models and/or couplings.

## Creating a New Model

In this section, we describe in detail how to construct a new model.  To illustrate this process, we also provide example code which could be used to re-implement the [`TypicalSection`](@ref) model.

### Manipulating a Model's Governing Equations

Before a model can be used with this package, its governing equations must be manipulated so that it satisfies the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x, y, p, t)`` is a function which defines the mass matrix corresponding to the differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of state variables, ``y`` is a vector of inputs, ``p`` is a vector of parameters, and ``t`` is the current time.  State variables are variables which have rate equations associated with them.  Inputs, or coupling variables, are variables  which may be defined in time using other models.  Parameters are variables which are user-specified and constant in time.

For example, the governing differential equations for the [`TypicalSection`](@ref) model are typically expressed as the second order ordinary differential equation
```math
\begin{bmatrix} m & S_\theta \\ S_\theta & I_\theta \end{bmatrix}
\begin{Bmatrix} \ddot{h} \\ \ddot{\theta} \end{Bmatrix} +
\begin{bmatrix} k_h & 0 \\ 0 & k_h \end{bmatrix}
\begin{Bmatrix} h \\ \theta \end{Bmatrix} =
\begin{Bmatrix} -\mathcal{L} \\ \mathcal{M} \end{Bmatrix}
```
where ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the structural imbalance, ``I_\theta`` is the mass moment of inertia, ``\mathcal{L}`` is the lift per unit span, and ``\mathcal{M}`` is the moment per unit span.  Expressed in the form expected by this package, the governing differential equation for the typical section model is
```math
M \dot{x} = K x + D y
```
where
```math
x = \begin{bmatrix} h & \theta & \dot{h} & \dot{\theta} \end{bmatrix}^T \quad
y = \begin{bmatrix} \mathcal{L} & \mathcal{M} \end{bmatrix}^T
```
```math
M =
\begin{bmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & m & S_\theta \\
    0 & 0 & S_\theta & I_P
\end{bmatrix}
\quad
K =
\begin{bmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
-k_h & 0 & 0 & 0 \\
0 & -k_\theta & 0 & 0
\end{bmatrix}
\quad
D =
\begin{bmatrix}
0 & 0 \\
0 & 0 \\
-1 & 0 \\
0 & 1
\end{bmatrix}
```

A special type of model which trivially satisfies the form of the governing differential equations expected by this package is a model with no state variables and/or inputs.  These models are designated as being subtypes of abstract type [`NoStateModel`](@ref) and are used solely to define the inputs of other models.  For example, the [`Steady`](@ref) and [`QuasiSteady`](@ref) models may be used to calculate the inputs corresponding to the [`TypicalSection`](@ref) model, but have no state variables of their own.

### Defining a Model's Type

Once the governing differential equations for a model has been manipulated into the expected format, and the state variables, inputs, and parameters for the model have been identified, the first step in defining a new model is to define a new type.  A docstring should also be provided along with the new type definition which defines the identities of the elements of the state, input, and parameter vectors for the model.  For example, the type definition for the [`TypicalSection`](@ref) model is
```julia
"""
    TypicalSection <: AbstractModel

Typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
struct TypicalSection <: AbstractModel end
```

### Defining Model Traits

The next step in defining a new model is to define the model's properties.  At a minimum, this requires defining new methods for the [`number_of_states`](@ref), [`number_of_inputs`](@ref), [`number_of_parameters`](@ref), and [`inplaceness`](@ref) functions, though additional method definitions may be necessary.  For example, the [`TypicalSection`](@ref) model's properties are defined in the following block of code.

```julia
number_of_states(::Type{TypicalSection}) = 4
number_of_inputs(::Type{TypicalSection}) = 2
number_of_parameters(::Type{TypicalSection}) = 5
inplaceness(::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{TypicalSection}) = Linear()
input_jacobian_type(::Type{TypicalSection}) = Constant()
```

The length of the model's state, input and parameter vectors is specified by defining new methods for the [`number_of_states`](@ref), [`number_of_inputs`](@ref), and [`number_of_parameters`](@ref) functions.  For out-of-place models, these methods must operate the model type so that the vector sizes are completely inferrable.  For in-place models, this restriction is loosened and these methods may operate on model instances instead.

Whether a model uses in-place or out-of-place function definitions is specified by defining a new method for the [`inplaceness`](@ref) function, which operates on the model type. For performance reasons, in-place functions are generally preferred.  The one exception is for models with small numbers of state variables, in which case the preferred approach is to use static arrays with out-of-place functions.

The properties of the model's mass matrix, state jacobian, and/or input jacobian are defined by defining new methods for the [`mass_matrix_type`](@ref), [`state_jacobian_type`](@ref), and/or [`input_jacobian_type`](@ref) functions, respectively.  By default, these properties assume their loosest possible definitions.

### Defining Methods for Governing Equations

Once the properties of a model have been defined, methods must be provided for the model which define its governing equations.  The right hand side of the governing structural differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models or the [`get_rates!`](@ref) function for in-place models.  For models with mass matrices, a new method must also be defined for the [`get_mass_matrix`](@ref) function (or [`get_mass_matrix!`](@ref) function if the model's functions are in-place functions).  For constant mass matrices (`mass_matrix_type(typeof(model)) == Constant()`), this function should be defined without the `x`, `y`, `p`, and `t` arguments.  

For example, the governing equations for the [`TypicalSection`](@ref) model may be defined using the following block of code

```julia
function get_rates(::TypicalSection, x, y, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract inputs
    L, M = y
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate state rates
    return SVector(hdot, θdot, -kh*h - L, -kθ*θ + M)
end

function get_mass_matrix(::MyTypicalSection, x, y, p, t)
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate mass matrix
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m Sθ; 0 0 Sθ Iθ]
end
```

### Defining Methods for Jacobians

Unless otherwise specified, the jacobian of the governing differential equations for a given model with respect to the state variables and/or inputs is calculated when necessary using forward automatic differentiation (enabled by the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package).  While this approach for computing the jacobians is convenient and exact, alternative methods for computing jacobians may be more computationally efficient.  To manually define the jacobian of the right hand side of the governing equations with respect to the state variables, a new method for [`get_state_jacobian`](@ref) (or [`get_state_jacobian!`](@ref) for in-place models) may be defined.  To manually define the jacobian of the right hand side of the governing equations with respect to the inputs, a new method for [`get_input_jacobian`](@ref) may be defined.  For example, for the [`TypicalSection`](@ref) model, these jacobians may be defined analytically using the following block of code

```julia
function get_state_jacobian(::MyTypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

get_input_jacobian(::MyTypicalSection) = @SMatrix [0 0; 0 0; -1 0; 0 1]
```

Note that while there is no in-place version of the [`get_input_jacobian`](@ref) function it may be constructed as either a linear map (for large matrices) or static array (for small matrices) in order to avoid allocations.

### Defining Methods for Unit Testing

In order to test whether provided mass matrices are correct for a given model, a new method for [`get_lhs()`](@ref) (which defines the left hand side of the governing differential equations) must be provided.  Since this function is only used for testing, there is no in-place version of this function.  For the [`TypicalSection`](@ref) model this function could be defined as follows

```julia
function get_lhs(::TypicalSection, dq, q, r, p, t)
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # extract state rates
    dh, dθ, dhdot, dθdot = dq
    # calculate mass matrix product
    return SVector(dh, dθ, m*dhdot + Sθ*dθdot, Sθ*dhdot + Iθ*dθdot)
end
```

## Creating a New Model Coupling

In this section, we describe in detail how to construct a new model by coupling multiple existing models together.  To illustrate this process, we also provide example code which could be used to re-implement the [`Wagner`](@ref) + [`TypicalSection`](@ref) model coupling.

### Coupled Model Theory

Coupled models in this package are defined by concatenating the governing differential equations, state variables, inputs, and parameters of multiple independent models into a single monolithic system.  As a result, coupled models defined by this package have the same general form as independent models in this package, that is, they may be defined by the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x, y, p, t)`` is a function which defines the mass matrix corresponding to the coupled model's governing differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of states, ``y`` is a vector of inputs (coupling variables), ``p`` is a vector of parameters, and ``t`` is the current time.

We assume the inputs for coupled models may be expressed as a function of the coupled model's state variables and parameters, as well as the current time.  Denoting such a function as `g(x, p, t)`, the inputs may then be expressed as ``y = g(x, p, t)``.  The governing differential equations for the coupled model then reduce to
```math
M(x, p, t)\dot{x} = f(x, p, t)
```


For example, the

This allows each model to influence other models through

If we now assume that the inputs to the coupled model may be defined as a function of the state variables and parameters of the coupled model, as well as the current time, then

We also assume the inputs for the coupled model may be expressed as a function of the coupled model's state variables and parameters, as well as the current time.  This assumption allows the state variables and parameters from each coupled model to influence the governing differential equations of the other coupled models through their input variables.

To model the interactions between models, we assume that the inputs to each component model may be defined as a function of the state variables, par



we assume that each model the interactions between models may be captured through each model's influence on the inputs corresponding to the governing equations of other models.

We also assume the inputs for the coupled model may be expressed as a function of the coupled model's state variables and parameters, as well as the current time.  This assumption allows the state variables and parameters from each coupled model to influence the governing differential equations of the other coupled models through their input variables.

We assume the inputs for the coupled model may be expressed


Each of the models which make up a coupled model influence other models through the other models inputs.  We


The inputs which correspond to each
Each of the models which make up

```math
y = g(x, p, t) - M_y(x,p,t) \dot{x} =
\begin{Bmatrix}
g_\text{aero}(x, p, t) - M_{y,\text{aero}}(x,p,t) \dot{x} \\
g_\text{stru}(x, p, t) - M_{y,\text{stru}}(x,p,t) \dot{x}
\end{Bmatrix}
```


To model the interactions between models

Coupled models defined by this package also satisfy the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x, y, p, t)`` is a function which defines the mass matrix corresponding to the differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of states, ``y`` is a vector of inputs (coupling variables), ``p`` is a vector of parameters, and ``t`` is the current time.

The state variables for coupled models

## Aeroelastic Model

To form a new aeroelastic model, we will be coupling Peters' finite state model with the typical section structural model.

### Theory

State variables for all models (including coupled models) in this package satisfy the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x, y, p, t)`` is a function which defines the mass matrix corresponding to the differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of states, ``y`` is a vector of inputs (coupling variables), ``p`` is a vector of parameters, and ``t`` is the current time.

We assume the state variables of the coupled model correspond to the state variables of all of its subcomponent models concatenated, the inputs of the coupled model correspond to the inputs of all of its subcomponent models concatenated, and the parameters of the coupled model correspond to the parameters of all of its subcomponent models concatenated, followed by an additional set of parameters.

```math
x = \begin{Bmatrix} x_\text{aero} \\ x_\text{stru} \end{Bmatrix} \quad
y = \begin{Bmatrix} y_\text{aero} \\ y_\text{stru} \end{Bmatrix} \quad
p = \begin{Bmatrix} p_\text{aero} \\ p_\text{stru} \\ p_\text{additional} \end{Bmatrix}
```

We also assume the inputs for the coupled model may be expressed as a function of the coupled model's state variables and parameters, as well as the current time.  This assumption allows the state variables and parameters from each coupled model to influence the governing differential equations of the other coupled models through their input variables.

```math
y = g(x, p, t) = \begin{Bmatrix} g_\text{aero}(x, p, t) \\ g_\text{stru}(x, p, t) \end{Bmatrix}
```

If the state rates of the coupled model are linearly dependent on the inputs, the definition of the inputs of the coupled model may be expanded to include a linear dependency on the coupled model's state rates.  Since this term will eventually be moved to the left side of the governing differential equations for the coupled model, we define it as a negative quantity.

```math
y = g(x, p, t) - M_y(x,p,t) \dot{x} =
\begin{Bmatrix}
g_\text{aero}(x, p, t) - M_{y,\text{aero}}(x,p,t) \dot{x} \\
g_\text{stru}(x, p, t) - M_{y,\text{stru}}(x,p,t) \dot{x}
\end{Bmatrix}
```

With these assumptions, the coupled system of equations may be defined as
```math
M(x,p,t) \dot{x} = f(x,p,t)
```
```math
M = \begin{bmatrix}
    M_\text{aero} & 0 \\
    0 & M_\text{stru}
\end{bmatrix} +
\begin{bmatrix}
    \frac{\partial f_\text{aero}}{\partial y_\text{aero}} & 0 \\
    0 & \frac{\partial f_\text{stru}}{\partial y_\text{stru}}
\end{bmatrix}
\begin{bmatrix}
    M_{y,\text{aero}} \\
    M_{y,\text{stru}}
\end{bmatrix} \quad
```
and the associated jacobian is
```math
\frac{\partial f_x}{\partial x} = \begin{bmatrix}
\frac{\partial f_\text{aero}}{\partial x_\text{aero}} & 0 \\
0 & \frac{\partial f_\text{stru}}{\partial x_\text{stru}}
\end{bmatrix} +
\begin{bmatrix}
    \frac{\partial f_\text{aero}}{\partial y_\text{aero}} & 0 \\
    0 & \frac{\partial f_\text{stru}}{\partial y_\text{stru}}
\end{bmatrix}
\begin{bmatrix}
    \frac{\partial g_\text{aero}}{\partial x_\text{aero}} &
    \frac{\partial g_\text{aero}}{\partial x_\text{stru}} \\
    \frac{\partial g_\text{stru}}{\partial x_\text{aero}} &
    \frac{\partial g_\text{stru}}{\partial x_\text{stru}}
\end{bmatrix}
```

#### Aerodynamic Model Inputs

The freestream velocity components ``u`` and ``v`` are assumed to be aligned with the undeflected chordwise and normal directions, respectively, so that
```math
u \approx U_\infty \\
v \approx \dot{h} \\
\omega \approx \dot{\theta}
```
where ``U_\infty`` is the freestream velocity magnitude, ``\theta`` is pitch, and ``h`` is plunge. To capture the effect of twist on the circulatory lift (since it is no longer implicitly modeled by the ``\frac{v}{u}`` quantity) twist is added to the effective angle of attack from Peter's finite state model so that the effective angle of attack is now given by
```math
\alpha_\text{eff} = \theta - \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega  + \frac{\lambda_0}{u} - \alpha_0
```
The original expression for the effective angle of attack may be used by defining the new variable ``\bar{v} = u \theta + v`` such that
```math
\alpha_\text{eff} = -\frac{\bar{v}}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega + \frac{\lambda_0}{u} - \alpha_0
```

With these assumptions, the inputs for Peters' finite state model are
```math
u = U_\infty \\
\omega = \dot{\theta} \\
\dot{v} = \ddot{h} \\
\dot{\omega} = \ddot{\theta}
```

#### Structural Model Inputs

A small angle assumption is also used to define the lift about the reference location as
```math
\mathcal{L} \approx \mathcal{N}
```
where ``\mathcal{N}`` is the normal force per unit span at the reference location.

With this assumption, the inputs for the typical section model is
```math
\mathcal{L} = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u\omega - a b \dot{\omega} \right) \\
\mathcal{M} = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\omega + b \left( \frac{1}{8} - \frac{a}{2} \right) \dot{\omega} \right] + b \left(\frac{1}{2} + a \right) \mathcal{L}
```
where
```math
\alpha = \theta - \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega + \frac{\lambda_0}{u} - \alpha_0
```
and
```math
\lambda_0 \approx \frac{1}{2} \sum_{n=1}^N b_n \lambda_n
```

### Defining the Coupled Model

To define the coupled model, we define a new method for the [`couple_models`](@ref) function.  This method should just return a tuple of models in the order in which their state variables, inputs, and parameters are defined.

```@example developer
"""
    couple_models(aero::Peters, stru::TypicalSection)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al. and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::MyPeters, stru::MyTypicalSection) = (aero, stru)

nothing #hide
```

### Defining Input Function Properties

Before we define the combined input function (and its associated functions) we need to define a few of its properties.

To define whether the combined input function should use an in-place or out-of-place format, we define a new method for the [`inplaceness`](@ref) function.

```@example developer
inplaceness(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = OutOfPlace()
nothing #hide
```

To indicate input function mass matrix properties, we define a new method for the [`mass_matrix_type`](@ref) function.

```@example developer
mass_matrix_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Linear()
nothing #hide
```

To indicate the properties of the input function jacobian with respect to the state variables, we define a new method for the [`state_jacobian_type`](@ref) function.

```@example developer
state_jacobian_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Nonlinear()
nothing #hide
```

To define the number of additional parameters introduced by this model, we define a new method for the [`number_of_parameters`](@ref) function.

```@example developer
number_of_parameters(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = 2
nothing
```

### Defining Input Function Methods

The portion of the inputs which is independent of the state rates is calculated using the [`get_inputs`](@ref) function for out-of-place aerodynamic load calculations and [`get_inputs!`](@ref) function for in-place aerodynamic load calculations.

```@example developer

function get_inputs(aero::MyPeters{N,TF,SV,SA}, stru::MyTypicalSection,
    s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    h, θ, hdot, θdot = s[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = U
    vbar = U*θ + hdot
    ω = θdot
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(vbar + d*ω - λ0 - u*α0) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, 0, 0, L, M)
end
nothing #hide
```

For out-of-place input functions, the mass matrix is calculated using the [`AerostructuralDynamics.get_input_mass_matrix`](@ref) function.  For in-place combined input functions, the mass matrix is calculated using the [`AerostructuralDynamics.get_input_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `u`, `p`, and `t` arguments.

```@example developer

function get_input_mass_matrix(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    vdot_dhdot = 1
    ωdot_dθdot = 1
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_dhdot = tmp/b
    L_dθdot = -tmp*a
    # moment at reference point
    M_dhdot = -tmp/2 + (b/2 + a*b)*L_dhdot
    M_dθdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_dθdot
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 -vdot_dhdot 0; 0 0 0 -ωdot_dθdot]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = @SMatrix [0 0 -L_dhdot -L_dθdot; 0 0 -M_dhdot -M_dθdot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

nothing #hide
```

### Input Function Performance Overloads

The jacobian of ``f_y`` with respect to the state variables is calculated using the [`AerostructuralDynamics.get_input_state_jacobian`](@ref) function for out-of-place combined input functions and [`AerostructuralDynamics.get_input_state_jacobian!`](@ref) function for in-place combined inputs functions.

```@example developer
function get_input_state_jacobian(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    ω_θdot = 1
    # calculate aerodynamic loads
    L_λ = -a0*ρ*u*b/2*bbar'
    M_λ = (b/2 + a*b)*L_λ
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    L_hdot = a0*ρ*u*b
    M_hdot = (b/2 + a*b)*L_hdot
    L_θdot = a0*ρ*u*b*(b/2 - a*b) + pi*ρ*b^2*u
    M_θdot = -pi*ρ*b^3*u + (b/2 + a*b)*L_θdot
    # construct sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 0; 0 0 0 ω_θdot; 0 0 0 0; 0 0 0 0]
    Jra = [L_λ; M_λ]
    Jrs = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end
nothing #hide
```

### Coupled Model Code

Putting it all together, a complete representation of the input function associated with our coupled model may be defined using the following block of code.

```@example developer-combined
"""
    couple_models(aero::Peters, stru::TypicalSection)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al. and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::MyPeters, stru::MyTypicalSection) = (aero, stru)

# --- traits --- #

inplaceness(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Linear()
state_jacobian_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Nonlinear()
number_of_parameters(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = 2

# --- methods --- #

function get_inputs(aero::MyPeters{N,TF,SV,SA}, stru::MyTypicalSection,
    s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    h, θ, hdot, θdot = s[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = U
    vbar = U*θ + hdot
    ω = θdot
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(vbar + d*ω - λ0 - u*α0) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, 0, 0, L, M)
end


function get_input_mass_matrix(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    vdot_dhdot = 1
    ωdot_dθdot = 1
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_dhdot = tmp/b
    L_dθdot = -tmp*a
    # moment at reference point
    M_dhdot = -tmp/2 + (b/2 + a*b)*L_dhdot
    M_dθdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_dθdot
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 -vdot_dhdot 0; 0 0 0 -ωdot_dθdot]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = @SMatrix [0 0 -L_dhdot -L_dθdot; 0 0 -M_dhdot -M_dθdot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    ω_θdot = 1
    # calculate aerodynamic loads
    L_λ = -a0*ρ*u*b/2*bbar'
    M_λ = (b/2 + a*b)*L_λ
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    L_hdot = a0*ρ*u*b
    M_hdot = (b/2 + a*b)*L_hdot
    L_θdot = a0*ρ*u*b*(b/2 - a*b) + pi*ρ*b^2*u
    M_θdot = -pi*ρ*b^3*u + (b/2 + a*b)*L_θdot
    # construct sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 0; 0 0 0 ω_θdot; 0 0 0 0; 0 0 0 0]
    Jra = [L_λ; M_λ]
    Jrs = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end
```

### Example Usage

```@example developer
using AerostructuralDynamics, LinearAlgebra

# reduced velocity range
V = range(0, 3.1, length=5000) # = U/(b*ωθ) (reduced velocity)

# non-dimensional parameters
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle

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

# dimensionalized velocity
U = V*b*ωθ

# aerodynamic models
aerodynamic_model = MyPeters{6}()

# structural model
structural_model = MyTypicalSection()

# coupled model
model = couple_models(aerodynamic_model, structural_model)

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(model), length(V))

# loop through each reduced frequency
for i = 1:length(V)
    # state variables
    u_aero = zeros(number_of_states(aerodynamic_model))
    u_stru = zeros(number_of_states(structural_model))
    u = vcat(u_aero, u_stru)

    # parameters
    p_aero = [a, b, a0, α0]
    p_stru = [kh, kθ, m, Sθ, Iθ]
    p_input = [U[i], ρ]
    p = vcat(p_aero, p_stru, p_input)

    # time
    t = 0.0

    # calculate inputs
    y = get_inputs(model, u, p, t)

    # perform linear stability analysis
    λ[:,i], Uλ, Vλ = get_eigen(model, u, y, p, t)
end

nothing #hide
```

```@example developer
using Plots
pyplot()

default(
    titlefontsize = 14,
    legendfontsize = 11,
    guidefontsize = 14,
    tickfontsize = 11,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true,
    framestyle = :zerolines)

sp1 = plot(
    title = "Non-Dimensional Frequency",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    legend = :topright
    )

sp2 = plot(
    title = "Non-Dimensional Damping",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (-0.7, 0.605),
    ytick = -0.6:0.2:0.6,
    ylabel = "\$ \\frac{Γ}{\\omega_\\theta} \$",
    legend = :topleft
    )

for i = 1:size(λ, 1)
    scatter!(sp1, V, imag.(λ[i,:])/ωθ,
        label = "",
        color = 1,
        markersize = 1,
        markerstrokewidth = 0,
        )
end

for i = 1:size(λ, 1)
    scatter!(sp2, V, real.(λ[i,:])/ωθ,
        label = "",
        color = 1,
        markersize = 1,
        markerstrokewidth = 0,
        )
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800))

savefig(p1, "example-stability.svg") #hide

nothing #hide
```

![](example-stability.svg)

## Model Ordering

In general, we suggest that the following ordering of model state variables, inputs, and parameters is used when constructing input functions.

1. Aerodynamics
2. Structural
3. Rigid Body (when present)

## Avoiding Mass Matrices

In order to take advantage of as many features of the DifferentialEquations package as possible (including local sensitivity analysis), at this point in time we recommend that the governing differential equations for models be reformulated to avoid using mass matrices whenever possible.








In this guide, we define a new aerodynamic and structural model and then couple them together to form a new aeroelastic model.  The same process may be followed to define and/or couple any number of models together.

```@contents
Pages = ["developer.md"]
Depth = 3
```

```@setup developer
using AerostructuralDynamics, StaticArrays, LinearAlgebra
```

```@setup developer-combined
using AerostructuralDynamics, StaticArrays, LinearAlgebra
```

##

## Typical Section Model

The first model we will be constructing is the typical section structural model.  

![](typical-section.svg)

The equations of motion for this model are:
```math
\begin{bmatrix} m & S_\theta \\ S_\theta & I_\theta \end{bmatrix}
\begin{Bmatrix} \ddot{h} \\ \ddot{\theta} \end{Bmatrix} +
\begin{bmatrix} k_h & 0 \\ 0 & k_h \end{bmatrix}
\begin{Bmatrix} h \\ \theta \end{Bmatrix} =
\begin{Bmatrix} -\mathcal{L} \\ \mathcal{M} \end{Bmatrix}
```
where ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the structural imbalance, ``I_\theta`` is the mass moment of inertia, ``\mathcal{L}`` is the lift per unit span, and ``\mathcal{M}`` is the moment per unit span.

### Theory

State variables for all models in this package satisfy the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x, y, p, t)`` is a function which defines the mass matrix corresponding to the differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of states, ``y`` is a vector of inputs (coupling variables), ``p`` is a vector of parameters, and ``t`` is the current time.

The equations of motion for the typical section model when expressed in the form expected by this package are
```math
M \dot{x} = K x + D y
```
where
```math
x = \begin{bmatrix} h & \theta & \dot{h} & \dot{\theta} \end{bmatrix}^T \quad
y = \begin{bmatrix} -\mathcal{L} & \mathcal{M} \end{bmatrix}^T
```
```math
M =
\begin{bmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & m & m b x_\theta \\
    0 & 0 & m b x_\theta & I_P
\end{bmatrix}
\quad
K =
\begin{bmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
-k_h & 0 & 0 & 0 \\
0 & -k_\theta & 0 & 0
\end{bmatrix}
\quad
D =
\begin{bmatrix}
0 & 0 \\
0 & 0 \\
-1 & 0 \\
0 & 1
\end{bmatrix}
```

### Defining a New Type

We start creating our model by defining a new type.  

```@example developer
"""
    MyTypicalSection <: AbstractModel

Typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
struct MyTypicalSection <: AbstractModel end

nothing #hide
```

Note that the state, input, and parameter identities for each model should be documented in the docstring associated with their type since this docstring provides the primary source of documentation for a given model.

### Defining Model Properties

We now need to define a few model properties.

To indicate the number of state variables in our model, we define a new method for the [`number_of_states`](@ref) function.  In our case, we have four state variables ``h, \theta, \dot{h}, \dot{\theta}``.

```@example developer
AerostructuralDynamics.number_of_states(::Type{MyTypicalSection}) = 4
nothing #hide
```

To indicate the number of inputs (coupling variables), we define a new method for the [`number_of_inputs`](@ref) function.  In our case, we have two inputs ``\\mathcal{L}, \\mathcal{M}``.  We define these variables as inputs rather than parameters since they will be calculated as a function of the state variables of the coupled system, rather than prescribed directly.  For aeroelastic systems, the inputs for aerodynamic models typically correspond to structural deflections while the inputs for structural models typically correspond to aerodynamic loads.

```@example developer
AerostructuralDynamics.number_of_inputs(::Type{MyTypicalSection}) = 2
nothing #hide
```

To indicate the number of parameters, we define a new method for the [`number_of_parameters`](@ref) function.  In our case, we have five parameters ``k_h, k_\\theta, m, S_\\theta, I_\\theta``.

```@example developer
AerostructuralDynamics.number_of_parameters(::Type{MyTypicalSection}) = 5
nothing #hide
```

To indicate whether our model uses in-place or out-of-place functions, we define a new method for the [`inplaceness`](@ref) function.  For performance reasons, in-place functions are generally preferred.  The one exception is for models with small numbers of state variables, in which case the preferred approach is to use static arrays with out-of-place functions.  For this model, we use the latter approach.

```@example developer
inplaceness(::Type{MyTypicalSection}) = OutOfPlace()
nothing #hide
```

To indicate mass matrix properties, we define a new method for the [`mass_matrix_type`](@ref) function.

```@example developer
mass_matrix_type(::Type{MyTypicalSection}) = Linear()
nothing #hide
```

To indicate the properties of the jacobian of the mass matrix multiplied state rates with respect to the state variables, we define a new method for the [`state_jacobian_type`](@ref) function.  This method definition is only required if the jacobian with respect to the state variables is manually defined.

```@example developer
state_jacobian_type(::Type{MyTypicalSection}) = Linear()
nothing #hide
```

To indicate the properties of the jacobian of the mass matrix multiplied state rates with respect to the inputs, we define a new method for the [`input_jacobian_type`](@ref) function.  This method definition is only required if the jacobian with respect to the inputs is manually defined.

```@example developer
input_jacobian_type(::Type{MyTypicalSection}) = Constant()
nothing #hide
```

### Defining Model Methods

Now that we have defined the properties of our model, we need to define the governing equations for its state variables.

The right hand side of the governing structural differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models or the [`get_rates!`](@ref) function for in-place models.  

```@example developer
function get_rates(::MyTypicalSection, x, y, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract inputs
    L, M = y
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate state rates
    return SVector(hdot, θdot, -kh*h - L, -kθ*θ + M)
end

nothing #hide
```

Since our model uses a mass matrix, we also need to define a new method for the [`get_mass_matrix`](@ref) function (or [`get_mass_matrix!`](@ref) function if the model's functions are in-place functions).  For constant mass matrices, these functions are called without the `x`, `y`, `p`, and `t` arguments.  

```@example developer
function get_mass_matrix(::MyTypicalSection, x, y, p, t)
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate mass matrix
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m Sθ; 0 0 Sθ Iθ]
end
nothing #hide
```

## Performance Overloads

The code we have presented so far fully defines the governing equations for the structural state variables of the typical section model.  At this point, we can either provide custom definitions for the jacobians associated with these governing equations, or allow the jacobians to be calculated using forward automatic differentiation (enabled by the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package).  While the latter approach is convenient, the former approach is typically less computationally expensive.

The jacobian of the right hand side of the governing equations with respect to the state variables may be defined using the [`get_state_jacobian`](@ref) function for out-of-place models or [`get_state_jacobian!`](@ref) function for in-place models.

```@example developer
function get_state_jacobian(::MyTypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end
nothing #hide
```

The jacobian of the right hand side of the governing equations with respect to the inputs may be defined using the [`get_input_jacobian`](@ref) function.  There is no out-of-place form for this function, however, it may be constructed as either a linear map (for large matrices) or static array (for small matrices) in order to avoid allocations.

```@example developer
get_input_jacobian(::MyTypicalSection) = @SMatrix [0 0; 0 0; -1 0; 0 1]
nothing #hide
```

### Typical Section Model Code

Putting it all together, a complete representation of our typical section model for use with this package may be defined using the following block of code.

```@example developer-combined
"""
    MyTypicalSection <: AbstractModel

Typical section structural model with state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``, inputs ``\\mathcal{L}, \\mathcal{M}``, and parameters ``k_h,
k_\\theta, m, S_\\theta, I_\\theta``
"""
struct MyTypicalSection <: AbstractModel end

# --- Traits --- #

number_of_states(::Type{MyTypicalSection}) = 4
number_of_inputs(::Type{MyTypicalSection}) = 2
number_of_parameters(::Type{MyTypicalSection}) = 5
inplaceness(::Type{MyTypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{MyTypicalSection}) = Linear()
state_jacobian_type(::Type{MyTypicalSection}) = Linear()
input_jacobian_type(::Type{MyTypicalSection}) = Constant()

# --- Methods --- #

function get_rates(::MyTypicalSection, x, y, p, t)
    # extract state variables
    h, θ, hdot, θdot = x
    # extract inputs
    L, M = y
    # extract parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate state rates
    return SVector(hdot, θdot, -kh*h - L, -kθ*θ + M)
end

function get_mass_matrix(::MyTypicalSection, x, y, p, t)
    # extract structural parameters
    kh, kθ, m, Sθ, Iθ = p
    # calculate mass matrix
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m Sθ; 0 0 Sθ Iθ]
end

# --- Performance Overloads --- #

function get_state_jacobian(::MyTypicalSection, x, y, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

get_input_jacobian(::MyTypicalSection) = @SMatrix [0 0; 0 0; -1 0; 0 1]

nothing #hide
```

## Peters' Finite State Model

The second model we will be constructing is Peters' finite state unsteady aerodynamics model.

![](airfoil.svg)

The normal force and moment at the reference location, as calculated using Peters' finite state model are
```math
\mathcal{N} = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u\dot{\theta} - a b \ddot{\theta} \right) \\
\mathcal{M} = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right] + b \left(\frac{1}{2} + a \right) \mathcal{N}
```
where ``u`` is the local freestream velocity in the chordwise direction, ``v`` is the local freestream velocity in the normal direction, ``\omega`` is the freestream angular velocity, ``a`` defines the reference location, ``b`` is the semichord, ``a_0`` is the lift curve slope, ``\rho_\infty`` is the air density, and ``\alpha_\text{eff}`` is the effective angle of attack.  The effective angle of attack for this model is defined as
```math
\alpha_\text{eff} = -\frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \dot{\theta} + \frac{\lambda_0}{u}  - \alpha_0
```
where ``\lambda_0`` is the induced velocity, and ``\alpha_0`` is the zero lift angle of attack.

The induced velocity ``\lambda_0`` is approximated from a set of N induced-flow states ``\lambda_1, \lambda_2, \dots, \lambda_N`` as
```math
\lambda \approx \frac{1}{2} \sum_{n=1}^N b_n \lambda_n
```
The set of N first-order ordinary differential equations which govern the N finite aerodynamic states are derived by Peters et al. as
```math
\bar{A} \lambda + \frac{u}{b} \lambda = \bar{c} \left[ -\dot{v} + u\dot{\theta} + b \left(\frac{1}{2} - a \right) \ddot{\theta} \right]
```
where
```math
\bar{A} = \bar{D} + \bar{d} \bar{b}^T + \bar{c} \bar{d}^T + \frac{1}{2} \bar{c}  \bar{b}^T \\
\bar{D}_{nm} = \begin{cases}
\frac{1}{2n} & n=m+1 \\
\frac{-1}{2n} & n=m-1 \\
0 & n \neq m \pm 1 \\
\end{cases}
\quad
\bar{b}_n = \begin{cases}
\left( -1 \right)^{n-1} \frac{(N+n-1)!}{(N-n-1)!}\frac{1}{\left(n!\right)^2} & n \neq N \\
\left( -1 \right)^{n-1} & n = N
\end{cases}
\quad
\bar{c}_n = \frac{2}{n}
\quad
\bar{d}_n = \begin{cases}
\frac{1}{2} & n = 1 \\
0 & n \neq 1
\end{cases}
```

The equations of motion for Peters' finite state model, when expressed in the form expected by this package are:
```math
\bar{A}\dot{\lambda} = \bar{c}\left[ \dot{v} + u\omega + b \left(\frac{1}{2} - a\right) \dot{\omega} \right] - \frac{u}{b}\lambda
```

### Defining a New Type

We start creating our aerodynamic model by defining a new type.  

```@example developer
"""
    MyPeters{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``u, \\omega,
\\dot{v}, \\dot{\\omega}`` and parameters ``a, b, a_0, \\alpha_0``
"""
struct MyPeters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

nothing #hide
```

Here `N` is the number of aerodynamic state variables and `TF` is the floating point type used to represent the constant matrices/vectors ``\bar{A}``, ``\bar{b}``, and ``\bar{c}``.  We store these matrices/vectors as fields of the newly defined type since they are constants.

For convenience, we create a constructor which initializes matrix ``\bar{A}`` and vectors ``\bar{b}`` and ``\bar{c}`` given the number of aerodynamic state variables and floating point type.

```@example developer
"""
    MyPeters{N,TF=Float64}()

Initialize an object of type `MyPeters` which has `N` aerodynamic
degrees of freedom.
"""
MyPeters{N}() where N = MyPeters{N,Float64}()

function MyPeters{N,TF}() where {N,TF}

    b = zeros(TF, N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(big(N + n - 1))/factorial(big(N - n - 1))*
            1/factorial(big(n))^2
    end
    b[N] = (-1)^(N-1)

    c = zeros(TF, N)
    for n = 1:N
        c[n] = 2/n
    end

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[n, m] = 1/(2*n)
    end
    for m in 2:N
        n = m - 1
        D[n, m] = -1/(2*n)
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return MyPeters(SMatrix{N,N,TF}(A), SVector{N,TF}(b), SVector{N,TF}(c))
end

nothing #hide
```

### Defining Model Properties

We now need to define a few model properties.

To indicate the number of state variables in our model, we define a new method for the [`number_of_states`](@ref) function.  In this case, we have an arbitrary number of aerodynamic states, though typically 3-10 aerodynamic states are used.  This method must be defined for all models.

```@example developer
number_of_states(::Type{MyPeters{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
nothing #hide
```

To indicate the number of inputs, we define a new method for the [`number_of_inputs`](@ref) function.  In this case, we have four inputs ``u, \dot{\theta}, \ddot{h}, \ddot{\theta}``.  This method must be defined for all models.

```@example developer
number_of_inputs(::Type{<:MyPeters}) = 4
nothing #hide
```

To indicate the number of parameters, we define a new method for the [`number_of_parameters`](@ref) function.  In this case, we have four parameters ``a, b, a_0, u_0``.  This method must be defined for all models.

```@example developer
number_of_parameters(::Type{<:MyPeters}) = 4
nothing #hide
```

To indicate whether we plan to use in-place or out-of-place functions, we define a new method for the [`inplaceness`](@ref) function.  This method must be defined for all models.

```@example developer
inplaceness(::Type{<:MyPeters}) = OutOfPlace()
nothing #hide
```

To indicate mass matrix properties, we define a new method for the [`mass_matrix_type`](@ref) function.  This method must be defined for all models.

```@example developer
mass_matrix_type(::Type{<:MyPeters}) = Constant()
nothing #hide
```

To indicate the properties of the jacobian of the state rates with respect to the state variables, we define a new method for the [`state_jacobian_type`](@ref) function.

```@example developer
state_jacobian_type(::Type{<:MyPeters}) = Linear()
nothing #hide
```

To indicate the properties of the jacobian of the state rates with respect to the inputs, we define a new method for the [`input_jacobian_type`](@ref) function.

```@example developer
input_jacobian_type(::Type{<:MyPeters}) = Nonlinear()
nothing #hide
```

### Defining Model Methods

Now that we have defined the properties of our model, we need to define the governing equations for its state variables.

The right hand side of the governing differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models and [`get_rates!`](@ref) function for in-place models.  

```@example developer
function get_rates(model::MyPeters{N,TF,SV,SA}, x, y, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return cbar*(vdot + u*ω + (b/2-a*b)*ωdot) - u/b*λ
end

nothing #hide
```

Since our model uses a mass matrix, we also need to define a new method for the [`get_mass_matrix`](@ref) function (or [`get_mass_matrix!`](@ref) if the model's functions are in-place functions).  For constant mass matrices, this function is defined without the `x`, `y`, `p`, and `t` arguments.  

```@example developer
get_mass_matrix(model::MyPeters) = model.A
nothing #hide
```

### Performance Overloads

The code we have presented so far fully defines the governing equations for the aerodynamic state variables of the Peters' finite state model.  At this point, we can either provide custom definitions for the jacobians associated with these governing equations, or allow the jacobians to be calculated using forward automatic differentiation (enabled by the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package).

The jacobian of the right hand side of the governing equations with respect to the state variables may be defined using the [`get_state_jacobian`](@ref) function for out-of-place models and [`get_state_jacobian!`](@ref) function for in-place models.

```@example developer
function get_state_jacobian(model::MyPeters, x, y, p, t)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return -u/b*Diagonal(one.(cbar))
end

nothing #hide
```

The jacobian of the right hand side of the governing equations with respect to the inputs is defined using the [`get_input_jacobian`](@ref) function.  

```@example developer
function get_input_jacobian(model::MyPeters, x, y, p, t)
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return hcat(cbar*ω - λ/b, u*cbar, cbar, (b/2-a*b)*cbar)
end

nothing #hide
```

### Aerodynamic Model Code

Putting it all together, a complete representation Peters' finite state aerodynamic model may be defined using the following block of code.

```@example developer-combined
"""
    MyPeters{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``u, \\omega,
\\dot{v}, \\dot{\\omega}`` and parameters ``a, b, a_0, \\alpha_0``
"""
struct MyPeters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

# --- Constructors --- #

"""
    MyPeters{N,TF=Float64}()

Initialize an object of type `MyPeters` which has `N` aerodynamic
degrees of freedom.
"""
MyPeters{N}() where N = MyPeters{N,Float64}()

function MyPeters{N,TF}() where {N,TF}

    b = zeros(TF, N)
    for n = 1:N-1
        b[n] = (-1)^(n-1)*factorial(big(N + n - 1))/factorial(big(N - n - 1))*
            1/factorial(big(n))^2
    end
    b[N] = (-1)^(N-1)

    c = zeros(TF, N)
    for n = 1:N
        c[n] = 2/n
    end

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[n, m] = 1/(2*n)
    end
    for m in 2:N
        n = m - 1
        D[n, m] = -1/(2*n)
    end

    A = D + d*b' + c*d' + 1/2*c*b'

    return MyPeters(SMatrix{N,N,TF}(A), SVector{N,TF}(b), SVector{N,TF}(c))
end

# --- Traits --- #

number_of_states(::Type{MyPeters{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
number_of_inputs(::Type{<:MyPeters}) = 4
number_of_parameters(::Type{<:MyPeters}) = 4
inplaceness(::Type{<:MyPeters}) = OutOfPlace()
mass_matrix_type(::Type{<:MyPeters}) = Constant()
state_jacobian_type(::Type{<:MyPeters}) = Linear()
input_jacobian_type(::Type{<:MyPeters}) = Nonlinear()

# --- Methods --- #

function get_rates(model::MyPeters{N,TF,SV,SA}, x, y, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return cbar*(vdot + u*ω + (b/2-a*b)*ωdot) - u/b*λ
end

get_mass_matrix(model::MyPeters) = model.A

# --- Performance Overloads --- #

function get_state_jacobian(model::MyPeters, x, y, p, t)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return -u/b*Diagonal(one.(cbar))
end

function get_input_jacobian(model::MyPeters, x, y, p, t)
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(x)
    # extract inputs
    u, ω, vdot, ωdot = y
    # extract parameters
    a, b, a0, α0 = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return hcat(cbar*ω - λ/b, u*cbar, cbar, (b/2-a*b)*cbar)
end

nothing #hide
```

## Aeroelastic Model

To form a new aeroelastic model, we will be coupling Peters' finite state model with the typical section structural model.

### Theory

State variables for all models (including coupled models) in this package satisfy the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where ``M(x, y, p, t)`` is a function which defines the mass matrix corresponding to the differential equation, ``f(x, y, p, t)`` is a function which defines the mass matrix multiplied state rates, ``x`` is a vector of states, ``y`` is a vector of inputs (coupling variables), ``p`` is a vector of parameters, and ``t`` is the current time.

We assume the state variables of the coupled model correspond to the state variables of all of its subcomponent models concatenated, the inputs of the coupled model correspond to the inputs of all of its subcomponent models concatenated, and the parameters of the coupled model correspond to the parameters of all of its subcomponent models concatenated, followed by an additional set of parameters.

```math
x = \begin{Bmatrix} x_\text{aero} \\ x_\text{stru} \end{Bmatrix} \quad
y = \begin{Bmatrix} y_\text{aero} \\ y_\text{stru} \end{Bmatrix} \quad
p = \begin{Bmatrix} p_\text{aero} \\ p_\text{stru} \\ p_\text{additional} \end{Bmatrix}
```

We also assume the inputs for the coupled model may be expressed as a function of the coupled model's state variables and parameters, as well as the current time.  This assumption allows the state variables and parameters from each coupled model to influence the governing differential equations of the other coupled models through their input variables.

```math
y = g(x, p, t) = \begin{Bmatrix} g_\text{aero}(x, p, t) \\ g_\text{stru}(x, p, t) \end{Bmatrix}
```

If the state rates of the coupled model are linearly dependent on the inputs, the definition of the inputs of the coupled model may be expanded to include a linear dependency on the coupled model's state rates.  Since this term will eventually be moved to the left side of the governing differential equations for the coupled model, we define it as a negative quantity.

```math
y = g(x, p, t) - M_y(x,p,t) \dot{x} =
\begin{Bmatrix}
g_\text{aero}(x, p, t) - M_{y,\text{aero}}(x,p,t) \dot{x} \\
g_\text{stru}(x, p, t) - M_{y,\text{stru}}(x,p,t) \dot{x}
\end{Bmatrix}
```

With these assumptions, the coupled system of equations may be defined as
```math
M(x,p,t) \dot{x} = f(x,p,t)
```
```math
M = \begin{bmatrix}
    M_\text{aero} & 0 \\
    0 & M_\text{stru}
\end{bmatrix} +
\begin{bmatrix}
    \frac{\partial f_\text{aero}}{\partial y_\text{aero}} & 0 \\
    0 & \frac{\partial f_\text{stru}}{\partial y_\text{stru}}
\end{bmatrix}
\begin{bmatrix}
    M_{y,\text{aero}} \\
    M_{y,\text{stru}}
\end{bmatrix} \quad
```
and the associated jacobian is
```math
\frac{\partial f_x}{\partial x} = \begin{bmatrix}
\frac{\partial f_\text{aero}}{\partial x_\text{aero}} & 0 \\
0 & \frac{\partial f_\text{stru}}{\partial x_\text{stru}}
\end{bmatrix} +
\begin{bmatrix}
    \frac{\partial f_\text{aero}}{\partial y_\text{aero}} & 0 \\
    0 & \frac{\partial f_\text{stru}}{\partial y_\text{stru}}
\end{bmatrix}
\begin{bmatrix}
    \frac{\partial g_\text{aero}}{\partial x_\text{aero}} &
    \frac{\partial g_\text{aero}}{\partial x_\text{stru}} \\
    \frac{\partial g_\text{stru}}{\partial x_\text{aero}} &
    \frac{\partial g_\text{stru}}{\partial x_\text{stru}}
\end{bmatrix}
```

#### Aerodynamic Model Inputs

The freestream velocity components ``u`` and ``v`` are assumed to be aligned with the undeflected chordwise and normal directions, respectively, so that
```math
u \approx U_\infty \\
v \approx \dot{h} \\
\omega \approx \dot{\theta}
```
where ``U_\infty`` is the freestream velocity magnitude, ``\theta`` is pitch, and ``h`` is plunge. To capture the effect of twist on the circulatory lift (since it is no longer implicitly modeled by the ``\frac{v}{u}`` quantity) twist is added to the effective angle of attack from Peter's finite state model so that the effective angle of attack is now given by
```math
\alpha_\text{eff} = \theta - \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega  + \frac{\lambda_0}{u} - \alpha_0
```
The original expression for the effective angle of attack may be used by defining the new variable ``\bar{v} = u \theta + v`` such that
```math
\alpha_\text{eff} = -\frac{\bar{v}}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega + \frac{\lambda_0}{u} - \alpha_0
```

With these assumptions, the inputs for Peters' finite state model are
```math
u = U_\infty \\
\omega = \dot{\theta} \\
\dot{v} = \ddot{h} \\
\dot{\omega} = \ddot{\theta}
```

#### Structural Model Inputs

A small angle assumption is also used to define the lift about the reference location as
```math
\mathcal{L} \approx \mathcal{N}
```
where ``\mathcal{N}`` is the normal force per unit span at the reference location.

With this assumption, the inputs for the typical section model is
```math
\mathcal{L} = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u\omega - a b \dot{\omega} \right) \\
\mathcal{M} = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\omega + b \left( \frac{1}{8} - \frac{a}{2} \right) \dot{\omega} \right] + b \left(\frac{1}{2} + a \right) \mathcal{L}
```
where
```math
\alpha = \theta - \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega + \frac{\lambda_0}{u} - \alpha_0
```
and
```math
\lambda_0 \approx \frac{1}{2} \sum_{n=1}^N b_n \lambda_n
```

### Defining the Coupled Model

To define the coupled model, we define a new method for the [`couple_models`](@ref) function.  This method should just return a tuple of models in the order in which their state variables, inputs, and parameters are defined.

```@example developer
"""
    couple_models(aero::Peters, stru::TypicalSection)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al. and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::MyPeters, stru::MyTypicalSection) = (aero, stru)

nothing #hide
```

### Defining Input Function Properties

Before we define the combined input function (and its associated functions) we need to define a few of its properties.

To define whether the combined input function should use an in-place or out-of-place format, we define a new method for the [`inplaceness`](@ref) function.

```@example developer
inplaceness(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = OutOfPlace()
nothing #hide
```

To indicate input function mass matrix properties, we define a new method for the [`mass_matrix_type`](@ref) function.

```@example developer
mass_matrix_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Linear()
nothing #hide
```

To indicate the properties of the input function jacobian with respect to the state variables, we define a new method for the [`state_jacobian_type`](@ref) function.

```@example developer
state_jacobian_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Nonlinear()
nothing #hide
```

To define the number of additional parameters introduced by this model, we define a new method for the [`number_of_parameters`](@ref) function.

```@example developer
number_of_parameters(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = 2
nothing
```

### Defining Input Function Methods

The portion of the inputs which is independent of the state rates is calculated using the [`get_inputs`](@ref) function for out-of-place aerodynamic load calculations and [`get_inputs!`](@ref) function for in-place aerodynamic load calculations.

```@example developer

function get_inputs(aero::MyPeters{N,TF,SV,SA}, stru::MyTypicalSection,
    s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    h, θ, hdot, θdot = s[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = U
    vbar = U*θ + hdot
    ω = θdot
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(vbar + d*ω - λ0 - u*α0) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, 0, 0, L, M)
end
nothing #hide
```

For out-of-place input functions, the mass matrix is calculated using the [`AerostructuralDynamics.get_input_mass_matrix`](@ref) function.  For in-place combined input functions, the mass matrix is calculated using the [`AerostructuralDynamics.get_input_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `u`, `p`, and `t` arguments.

```@example developer

function get_input_mass_matrix(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    vdot_dhdot = 1
    ωdot_dθdot = 1
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_dhdot = tmp/b
    L_dθdot = -tmp*a
    # moment at reference point
    M_dhdot = -tmp/2 + (b/2 + a*b)*L_dhdot
    M_dθdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_dθdot
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 -vdot_dhdot 0; 0 0 0 -ωdot_dθdot]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = @SMatrix [0 0 -L_dhdot -L_dθdot; 0 0 -M_dhdot -M_dθdot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

nothing #hide
```

### Input Function Performance Overloads

The jacobian of ``f_y`` with respect to the state variables is calculated using the [`AerostructuralDynamics.get_input_state_jacobian`](@ref) function for out-of-place combined input functions and [`AerostructuralDynamics.get_input_state_jacobian!`](@ref) function for in-place combined inputs functions.

```@example developer
function get_input_state_jacobian(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    ω_θdot = 1
    # calculate aerodynamic loads
    L_λ = -a0*ρ*u*b/2*bbar'
    M_λ = (b/2 + a*b)*L_λ
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    L_hdot = a0*ρ*u*b
    M_hdot = (b/2 + a*b)*L_hdot
    L_θdot = a0*ρ*u*b*(b/2 - a*b) + pi*ρ*b^2*u
    M_θdot = -pi*ρ*b^3*u + (b/2 + a*b)*L_θdot
    # construct sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 0; 0 0 0 ω_θdot; 0 0 0 0; 0 0 0 0]
    Jra = [L_λ; M_λ]
    Jrs = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end
nothing #hide
```

### Coupled Model Code

Putting it all together, a complete representation of the input function associated with our coupled model may be defined using the following block of code.

```@example developer-combined
"""
    couple_models(aero::Peters, stru::TypicalSection)

Create an aerostructural model using the unsteady aerodynamic model defined by
Peters et al. and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::MyPeters, stru::MyTypicalSection) = (aero, stru)

# --- traits --- #

inplaceness(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Linear()
state_jacobian_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Nonlinear()
number_of_parameters(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = 2

# --- methods --- #

function get_inputs(aero::MyPeters{N,TF,SV,SA}, stru::MyTypicalSection,
    s, p, t) where {N,TF,SV,SA}
    # extract state variables
    λ = s[SVector{N}(1:N)]
    h, θ, hdot, θdot = s[SVector{4}(N+1:N+4)]
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    u = U
    vbar = U*θ + hdot
    ω = θdot
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # lift at reference point
    L = tmp1*(vbar + d*ω - λ0 - u*α0) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, ω, 0, 0, L, M)
end


function get_input_mass_matrix(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    vdot_dhdot = 1
    ωdot_dθdot = 1
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_dhdot = tmp/b
    L_dθdot = -tmp*a
    # moment at reference point
    M_dhdot = -tmp/2 + (b/2 + a*b)*L_dhdot
    M_dθdot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_dθdot
    # construct submatrices
    Mda = zeros(SMatrix{4,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 -vdot_dhdot 0; 0 0 0 -ωdot_dθdot]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = @SMatrix [0 0 -L_dhdot -L_dθdot; 0 0 -M_dhdot -M_dθdot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

# --- performance overloads --- #

function get_input_state_jacobian(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, s, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    bbar = aero.b
    # local freestream velocity components
    ω_θdot = 1
    # calculate aerodynamic loads
    L_λ = -a0*ρ*u*b/2*bbar'
    M_λ = (b/2 + a*b)*L_λ
    L_θ = a0*ρ*u^2*b
    M_θ = (b/2 + a*b)*L_θ
    L_hdot = a0*ρ*u*b
    M_hdot = (b/2 + a*b)*L_hdot
    L_θdot = a0*ρ*u*b*(b/2 - a*b) + pi*ρ*b^2*u
    M_θdot = -pi*ρ*b^3*u + (b/2 + a*b)*L_θdot
    # construct sub-matrices
    Jda = zeros(SMatrix{4,N,TF}) # d(d)/d(dλ)
    Jds = @SMatrix [0 0 0 0; 0 0 0 ω_θdot; 0 0 0 0; 0 0 0 0]
    Jra = [L_λ; M_λ]
    Jrs = @SMatrix [0 L_θ L_hdot L_θdot; 0 M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end
```

### Example Usage

```@example developer
using AerostructuralDynamics, LinearAlgebra

# reduced velocity range
V = range(0, 3.1, length=5000) # = U/(b*ωθ) (reduced velocity)

# non-dimensional parameters
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle

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

# dimensionalized velocity
U = V*b*ωθ

# aerodynamic models
aerodynamic_model = MyPeters{6}()

# structural model
structural_model = MyTypicalSection()

# coupled model
model = couple_models(aerodynamic_model, structural_model)

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(model), length(V))

# loop through each reduced frequency
for i = 1:length(V)
    # state variables
    u_aero = zeros(number_of_states(aerodynamic_model))
    u_stru = zeros(number_of_states(structural_model))
    u = vcat(u_aero, u_stru)

    # parameters
    p_aero = [a, b, a0, α0]
    p_stru = [kh, kθ, m, Sθ, Iθ]
    p_input = [U[i], ρ]
    p = vcat(p_aero, p_stru, p_input)

    # time
    t = 0.0

    # calculate inputs
    y = get_inputs(model, u, p, t)

    # perform linear stability analysis
    λ[:,i], Uλ, Vλ = get_eigen(model, u, y, p, t)
end

nothing #hide
```

```@example developer
using Plots
pyplot()

default(
    titlefontsize = 14,
    legendfontsize = 11,
    guidefontsize = 14,
    tickfontsize = 11,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true,
    framestyle = :zerolines)

sp1 = plot(
    title = "Non-Dimensional Frequency",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    legend = :topright
    )

sp2 = plot(
    title = "Non-Dimensional Damping",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (-0.7, 0.605),
    ytick = -0.6:0.2:0.6,
    ylabel = "\$ \\frac{Γ}{\\omega_\\theta} \$",
    legend = :topleft
    )

for i = 1:size(λ, 1)
    scatter!(sp1, V, imag.(λ[i,:])/ωθ,
        label = "",
        color = 1,
        markersize = 1,
        markerstrokewidth = 0,
        )
end

for i = 1:size(λ, 1)
    scatter!(sp2, V, real.(λ[i,:])/ωθ,
        label = "",
        color = 1,
        markersize = 1,
        markerstrokewidth = 0,
        )
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800))

savefig(p1, "example-stability.svg") #hide

nothing #hide
```

![](example-stability.svg)

## Model Ordering

In general, we suggest that the following ordering of model state variables, inputs, and parameters is used when constructing input functions.

1. Aerodynamics
2. Structural
3. Rigid Body (when present)

## Avoiding Mass Matrices

In order to take advantage of as many features of the DifferentialEquations package as possible (including local sensitivity analysis), at this point in time we recommend that the governing differential equations for models be reformulated to avoid using mass matrices whenever possible.
