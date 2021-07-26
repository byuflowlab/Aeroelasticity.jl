# Developer's Guide

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

Model parameters may be stored either as fields of the newly defined type or passed to the solver through the parameters vector.  In general, parameters which may be used as design variables as part of an optimization should be passed to the solver through the parameter vector, while constant parameters should be stored as fields of the newly defined type.  For this model, we define all parameters as elements of the parameter vector.

Note that the state, input, and parameter identities for each model should be documented in the docstring associated with their type since this docstring provides the primary source of documentation for a given model.

### Defining Model Properties

We now need to define a few model properties.

To indicate the number of state variables in our model, we define a new method for the [`number_of_states`](@ref) function.  In our case, we have four state variables ``h, \\theta, \\dot{h},
\\dot{\\theta}``.

```@example developer
number_of_states(::Type{MyTypicalSection}) = 4
nothing #hide
```

To indicate the number of inputs (coupling variables), we define a new method for the [`number_of_inputs`](@ref) function.  In our case, we have two inputs ``\\mathcal{L}, \\mathcal{M}``.  We define these variables as inputs rather than parameters since they will be calculated as a function of the state variables of the coupled system, rather than prescribed directly.  For aeroelastic systems, the inputs for aerodynamic models typically correspond to structural deflections while the inputs for structural models typically correspond to aerodynamic loads.

```@example developer
number_of_inputs(::Type{MyTypicalSection}) = 2
nothing #hide
```

To indicate the number of parameters, we define a new method for the [`number_of_parameters`](@ref) function.  In our case, we have five parameters ``k_h, k_\\theta, m, S_\\theta, I_\\theta``.

```@example developer
number_of_parameters(::Type{MyTypicalSection}) = 5
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

Since our model uses a mass matrix, we also need to define a new method for the [`get_mass_matrix`](@ref) function (or [`get_mass_matrix!`](@ref) function if the model's functions are in-place functions).  For constant mass matrices, these functions are called without the `λ`, `d`, `p`, and `t` arguments.  

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

The code we have presented so far fully defines the governing equations for the structural state variables of the typical section model.  At this point, we can either provide custom definitions for the jacobians associated with these governing equations, or allow the jacobians to be calculated using forward automatic differentiation enabled by the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package.  While the latter approach is convenient, the former approach typically allows for lower computational expenses when computing jacobians.

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

The jacobian of the right hand side of the governing equations with respect to the inputs may be defined using the [`get_input_jacobian`](@ref) function.  There is no out-of-place form for this function, however, it may be constructed as either a linear map (if large) or static array (if small) in order to avoid allocations.

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

Here `N` is the number of aerodynamic state variables and `TF` is the floating point type used to represent the matrices/vectors ``\bar{A}``, ``\bar{b}``, and ``\bar{c}``.  We store these matrices/vectors as fields of the newly defined type since they are constant for a given number of aerodynamic state variables.

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

To indicate the number of state variables in our model, we define a new method for the [`number_of_states`](@ref) function.  In our case, we have an arbitrary number of aerodynamic states ``\lambda``, though typically 3-10 aerodynamic states are used with Peters' finite state model.  This method must be defined for all models.

```@example developer
number_of_states(::Type{MyPeters{N,TF,SV,SA}}) where {N,TF,SV,SA} = N
nothing #hide
```

To indicate the number of inputs, we define a new method for the [`number_of_inputs`](@ref) function.  In our case, we have four inputs ``u, \dot{\theta}, \ddot{h}, \ddot{\theta}``.  This method must be defined for all models.

```@example developer
number_of_inputs(::Type{<:MyPeters}) = 4
nothing #hide
```

To indicate the number of parameters, we define a new method for the [`number_of_parameters`](@ref) function.  In our case, we have four parameters ``a, b, a_0, u_0``.  This method must be defined for all models.

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

The right hand side of the governing aerodynamic differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models and [`get_rates!`](@ref) function for in-place models.  

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

Since our model uses a mass matrix, we also need to define a new method for the [`get_mass_matrix`](@ref) function (or [`get_mass_matrix!`](@ref) if the model's functions are in-place functions).  For constant mass matrices, these functions are called without the `x`, `y`, `p`, and `t` arguments.  

```@example developer
get_mass_matrix(model::MyPeters) = model.A
nothing #hide
```

### Performance Overloads

The code we have presented so far fully defines the governing equations for the aerodynamic state variables of the Peters' finite state model.  At this point, we can either provide custom definitions for the jacobians associated with these governing equations, or allow the jacobians to be calculated using forward automatic differentiation enabled by the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package.  While the latter approach is convenient, the former approach typically allows for lower computational expenses when computing jacobians.

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

We assume the state variables of the coupled model correspond to the state variables of all of its subcomponent models concatenated, the inputs of the coupled model correspond to the inputs of all of its subcomponent models concatenated, and the parameters of the coupled model correspond to the parameters of all of its subcomponent models concatenated, followed by an additional set of parameters corresponding to the coupled model.

```math
x = \begin{bmatrix} x_\text{aero} & x_\text{stru} \end{bmatrix}^T \quad
y = \begin{bmatrix} y_\text{aero} & y_\text{stru} \end{bmatrix}^T \quad
p = \begin{bmatrix} p_\text{aero} & p_\text{stru} & p_\text{additional} \end{bmatrix}^T
```

We also assume the inputs for the coupled model may be expressed as a function of the coupled model's state variables and parameters, as well as the current time.  This assumption allows the state variables and parameters from each subcomponent model to influence the governing differential equations of the other subcomponent models through the input variables.

```math
y = g(x, p, t) = \begin{bmatrix} g_\text{aero}(x, p, t) \\ g_\text{stru}(x, p, t) \end{bmatrix}
```

If the state rates of the coupled model are linearly dependent on the inputs, the definition of the inputs of the coupled model may be expanded to include a linear dependency on the coupled model's state rates.  Since this term will eventually be moved to the left side of the governing differential equations for the coupled model, we define it with a negative sign.

```math
y = g(x, p, t) - M_y(x,p,t) \dot{x} =
\begin{bmatrix}
g_\text{aero}(x, p, t) - M_{y,\text{aero}}(x,p,t) \dot{x} \\
g_\text{stru}(x, p, t) - M_{y,\text{stru}}(x,p,t) \dot{x}
\end{bmatrix}
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

Introducing the additional parameter ``U_\infty``, representing the freestream velocity magnitude, we define the following variables for use with the aerodynamic model:

```math
u = U_\infty \\
v = U_\infty \theta + \dot{h}
\omega = \dot{theta}
\dot{v} = \ddot{h} \\
\dot{\omega} = \ddot{\theta}
```

Note that the quantity ``v`` is the vertical freestream velocity resulting from both the pitching angle and the unsteady plunging motion of the section.

#### Structural Model Inputs

The inputs to the structural model are ``
The aerodynamic loads calculated by Peters' finite state model when coupled with the typical section model may be expressed in the form expected by this package as
```math
r = J_{ra} \lambda + J_{rs} q - M_{rs} \dot{q}
```

where

```math
r =\begin{bmatrix}
L &
M_\frac{1}{4}
\end{bmatrix}^T
\quad
\lambda = \begin{bmatrix}
\lambda_1 &
\lambda_2 &
... &
\lambda_N
\end{bmatrix}^T
\quad
q = \begin{bmatrix}
h &
\theta &
\dot{h} &
\dot{\theta}
\end{bmatrix}^T
\\
J_{rs} = 2 \pi \rho_\infty b U
\begin{bmatrix}
0 & U & 1 &  \frac{b}{2} + b \left(\frac{1}{2} - a\right) \\
0 & 0 & 0 & - \frac{b^2}{2}
\end{bmatrix}
\quad
J_{ra} = - 2 \pi \rho_\infty b U \begin{bmatrix} \bar{b}^T \\ 0_{1 \times N} \end{bmatrix}
\\
M_{rs} = \pi \rho_\infty b^2
\begin{bmatrix}
0 & 0 & -1 & ba \\
0 & 0 & \frac{b}{2} & b^2 \left(\frac{1}{8} - \frac{a}{2}\right)
\end{bmatrix}
```

#### Combined Input Function

Upon combining the expressions for the aerodynamic loads and structural deflections, we obtain the following expression for the combined input function.

```math
y = J_{y} y - M_{y} \dot{u}
```
where
```math
J_{y} = \begin{bmatrix} 0 & J_{ds} \\ J_{ra} & J_{rs}  \end{bmatrix} \quad
M_{y} = \begin{bmatrix} 0 & M_{ds} \\ 0 & M_{rs} \end{bmatrix}
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
mass_matrix_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Varying()
nothing #hide
```

To indicate the properties of the input function jacobian with respect to the state variables, we define a new method for the [`state_jacobian_type`](@ref) function.

```@example developer
state_jacobian_type(::Type{<:MyPeters}, ::Type{MyTypicalSection}) = Varying()
nothing #hide
```

### Input Mass Matrix Equation

For out-of-place input functions, the mass matrix is calculated using the [`AerostructuralDynamics.get_input_mass_matrix`](@ref) function.  For in-place combined input functions, the mass matrix is calculated using the [`AerostructuralDynamics.get_input_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `u`, `p`, and `t` arguments.

```@example developer

function get_input_mass_matrix(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # construct submatrices
    Mda = zeros(SMatrix{3,N,TF})
    Mds = @SMatrix [0 0 0 0; 0 0 -1 0; 0 0 0 -1]
    Mra = zeros(SMatrix{2,N,TF})
    Mrs = hcat(
        zeros(SVector{2,TF}),
        zeros(SVector{2,TF}),
        -SVector(pi*ρ*b^2, -pi/2*ρ*b^3),
        -SVector(-pi*ρ*a*b^3, -pi/8*ρ*b^4*(1 - 4*a)))
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end

nothing #hide
```

### Input Equation

The portion of the inputs which is independent of the state rates is calculated using the [`get_inputs`](@ref) function for out-of-place aerodynamic load calculations and [`get_inputs!`](@ref) function for in-place aerodynamic load calculations.

```@example developer

function get_inputs(aero::MyPeters{N,TF,SV,SA}, stru::MyTypicalSection,
    u, p, t) where {N,TF,SV,SA}
    # indices for extracting state variables
    iλ = SVector{N}(1:N)
    iq = SVector{4}(N+1:N+4)
    # separate aerodynamic and structural states
    λ = u[iλ]
    q = u[iq]
    # extract structural state variables
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # extract model constants
    bbar = aero.b
    # calculate induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # calculate (partial) lift
    L = 2*pi*ρ*U*b*(hdot + U*θ + (b/2-a*b)*θdot - λ0) + pi*ρ*b^2*U*θdot
    # calculate (partial) quarter-chord moment
    M = -pi*ρ*b^3*U*θdot
    # return portion of inputs that is not dependent on the state rates
    return SVector(θdot, 0, 0, L, M)
end

nothing #hide
```

### Input Jacobian

The jacobian of ``f_y`` with respect to the state variables is calculated using the [`AerostructuralDynamics.get_input_state_jacobian`](@ref) function for out-of-place combined input functions and [`AerostructuralDynamics.get_input_state_jacobian!`](@ref) function for in-place combined inputs functions.

```@example developer
function get_input_state_jacobian(aero::MyPeters{N,TF,SV,SA},
    stru::MyTypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, U, ρ, a, b, kh, kθ, m, xθ, Ip = p
    # extract model constants
    bbar = aero.b
    # compute jacobian sub-matrices
    Jda = zeros(SMatrix{3,N,TF})
    Jds = @SMatrix [0 0 0 1; 0 0 0 0; 0 0 0 0]
    Jra = -pi*ρ*U*b*vcat(bbar', zero(bbar)')
    Jrs = hcat(
        SVector(0, 0),
        SVector(2*pi*ρ*U^2*b, 0),
        SVector(2*pi*ρ*U*b, 0),
        SVector(2*pi*ρ*b^2*U*(1 - a), -pi*ρ*b^3*U)
        )
    # return jacobian
    return [Jda Jds; Jra Jrs]
end
nothing #hide
```

## Model Ordering

In general, we suggest that the following ordering of model state variables, inputs, and parameters is used when constructing input functions.

1. Aerodynamics
2. Structural
3. Rigid Body (when present)

## Avoiding Mass Matrices

In order to take advantage of as many features of the DifferentialEquations package as possible (including local sensitivity analysis) we recommend that the governing differential equations for models be reformulated to avoid using mass matrices whenever possible.
