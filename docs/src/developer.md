# Developer's Guide

In this guide, we describe how to create new independent and/or coupled models.

```@contents
Pages = ["library.md"]
Depth = 3
```

## Theory

This package assumes that the governing equations for the state variables corresponding to all models satisfy the implicit ordinary differential equation
```math
0 = f(\dot{x},x,y,p,t)
```
where ``f(\dot{x}, x, y, p, t)`` is a residual function, ``x`` is a vector of state variables, ``y`` is a vector of inputs, ``p`` is a vector of parameters, and ``t`` is the current time.  In the context of this package, we define inputs as variables which may vary over time and parameters as variables which do not vary over time.  Typically, inputs are defined as any model parameter which could correspond to the output or outputs from other models.  We call t


The primary type of model encountered in this package is a model with state variables which are governed by  the following implicit ordinary differential equation.  

## Creating a New Model

In this section, we describe in detail how to construct a new independent model.  To demonstrate this process, we also show how one might implement the [`TypicalSection`](@ref) model.

### Manipulating a Model's Governing Equations

Before a model can be used with this package, its governing equations must be manipulated so that it satisfies the implicit ordinary differential equation
```math
0 = f(\dot{x},x,y,p,t)
```
where ``f(\dot{x}, x, y, p, t)`` is a residual function, ``x`` is a vector of state variables, ``y`` is a vector of inputs, ``p`` is a vector of parameters, and ``t`` is the current time.  State variables are variables which have rate equations associated with them.  Inputs are variables which may vary over time, possibly as defined by other models.  Parameters are variables which are user-specified and constant in time.

For example, the governing differential equation for the [`TypicalSection`](@ref) model is often expressed as the second order ordinary differential equation
```math
\begin{bmatrix} m & S_\theta \\ S_\theta & I_\theta \end{bmatrix}
\begin{bmatrix} \ddot{h} \\ \ddot{\theta} \end{bmatrix} +
\begin{bmatrix} k_h & 0 \\ 0 & k_h \end{bmatrix}
\begin{bmatrix} h \\ \theta \end{bmatrix} =
\begin{bmatrix} -\mathcal{L} \\ \mathcal{M} \end{bmatrix}
```
where ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``S_\theta`` is the structural imbalance, ``I_\theta`` is the mass moment of inertia, ``\mathcal{L}`` is the lift per unit span, and ``\mathcal{M}`` is the moment per unit span.  Expressed in the form expected by this package, the governing differential equation for this model is
```math
0 = M \dot{x} + K x + D y
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
0 & 0 & -1 & 0 \\
0 & 0 & 0 & -1 \\
k_h & 0 & 0 & 0 \\
0 & k_\theta & 0 & 0
\end{bmatrix}
\quad
D =
\begin{bmatrix}
0 & 0 \\
0 & 0 \\
1 & 0 \\
0 & -1
\end{bmatrix}
```

A special type of model which trivially satisfies the form of the governing differential equations expected by this package is a model with no state variables and/or inputs.  These models are designated as being subtypes of abstract type [`NoStateModel`](@ref AerostructuralDynamics.NoStateModel) and are used solely to define the inputs of other models.  For example, the [`Steady`](@ref) and [`QuasiSteady`](@ref) models may be used to calculate the inputs corresponding to the [`TypicalSection`](@ref) model, but have no state variables of their own.

### Defining a Model's Type

Once the governing differential equations for a model has been manipulated into the expected format, and the state variables, inputs, and parameters for the model have been identified, the first step in defining a new model is to define a new type.  A docstring should also be provided along with the new type definition which defines the identities of the elements of the state, input, and parameter vectors for the model.  

For example, the type definition for the [`TypicalSection`](@ref) is defined as
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

The next step in defining a new model is to define the model's properties.  At a minimum, this requires defining new methods for the [`number_of_states`](@ref), [`number_of_inputs`](@ref), [`number_of_parameters`](@ref), and [`inplaceness`](@ref AerostructuralDynamics.inplaceness) functions, though additional method definitions may be necessary.  

The length of the model's state, input and parameter vectors is specified by defining new methods for the [`number_of_states`](@ref), [`number_of_inputs`](@ref), and [`number_of_parameters`](@ref) functions.  For out-of-place models, these methods must operate the model type so that the vector sizes are completely inferrable.  For in-place models, this restriction is loosened and these methods may operate on model instances instead.

Whether a model uses in-place or out-of-place function definitions is specified by defining a new method for the [`inplaceness`](@ref AerostructuralDynamics.inplaceness) function, which operates on the model type. For performance reasons, in-place functions are generally preferred.  The one exception is for models with small numbers of state variables, in which case the preferred approach is to use static arrays with out-of-place functions.

The properties of the model's mass matrix, state jacobian, and/or input jacobian are defined by defining new methods for the [`mass_matrix_type`](@ref AerostructuralDynamics.mass_matrix_type), [`state_jacobian_type`](@ref AerostructuralDynamics.state_jacobian_type), and/or [`input_jacobian_type`](@ref AerostructuralDynamics.input_jacobian_type) functions, respectively.  By default, these properties assume their loosest possible definitions.

As an example, the [`TypicalSection`](@ref) model's properties are defined using the following block of code.

```julia
number_of_states(::Type{TypicalSection}) = 4
number_of_inputs(::Type{TypicalSection}) = 2
number_of_parameters(::Type{TypicalSection}) = 5
inplaceness(::Type{TypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{TypicalSection}) = Linear()
state_jacobian_type(::Type{TypicalSection}) = Linear()
input_jacobian_type(::Type{TypicalSection}) = Constant()
```

### Defining Methods for Governing Equations

Once the properties of a model have been defined, methods must be provided for the model which define its governing equations.  The right hand side of the governing structural differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models or the [`get_rates!`](@ref) function for in-place models.  For models with mass matrices, a new method must also be defined for the [`get_mass_matrix`](@ref) function (or [`get_mass_matrix!`](@ref) function if the model's functions are in-place).  For constant mass matrices (`mass_matrix_type(typeof(model)) == Constant()`), this function should be defined without the `x`, `y`, `p`, and `t` arguments.  

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

Unless otherwise specified, the jacobian of the governing differential equations for a given model with respect to the state variables and/or inputs is calculated when necessary using forward automatic differentiation (enabled by the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package).  While this approach for computing the jacobians is convenient and exact, alternative methods for computing jacobians may be more computationally efficient.  To manually define the jacobian of the right hand side of the governing equations with respect to the state variables, a new method for [`get_state_jacobian`](@ref) (or [`get_state_jacobian!`](@ref) for in-place models) may be defined.  To manually define the jacobian of the right hand side of the governing equations with respect to the inputs, a new method for [`get_input_jacobian`](@ref AerostructuralDynamics.get_input_jacobian) may be defined.  

For the [`TypicalSection`](@ref) model, these jacobians may be defined analytically using the following block of code

```julia
function get_state_jacobian(::MyTypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

get_input_jacobian(::MyTypicalSection) = @SMatrix [0 0; 0 0; -1 0; 0 1]
```

Note that while there is no in-place version of the [`get_input_jacobian`](@ref AerostructuralDynamics.get_input_jacobian) function it may be constructed as either a linear map (for large matrices) or static array (for small matrices) in order to avoid allocations.

### Defining Methods for Unit Testing

In order to test whether provided mass matrices are correct for a given model, a new method for [`get_lhs()`](@ref AerostructuralDynamics.get_lhs) (which defines the left hand side of the governing differential equations) must be provided.  Since this function is only used for testing, there is no in-place version of this function.  

For the [`TypicalSection`](@ref) model, the new method could be defined as follows

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

### Defining Convenience Methods

To aid users in defining state, input, and parameter vectors, new methods for the [`set_states!`](@ref), [`set_inputs!`](@ref), and [`set_parameters!`](@ref) functions should be provided.  For the [`TypicalSection`](@ref) model, these new methods may be defined using the following block of code.

```julia
function set_states!(x, model::TypicalSection; h, theta, hdot, thetadot)

    x[1] = h
    x[2] = theta
    x[3] = hdot
    x[4] = thetadot

    return x
end

function set_inputs!(y, model::TypicalSection; L, M)

    y[1] = L
    y[2] = M

    return y
end

function set_parameters!(p, model::TypicalSection; kh, ktheta, m, Stheta, Itheta)

    p[1] = kh
    p[2] = ktheta
    p[3] = m
    p[4] = Stheta
    p[5] = Itheta

    return p
end
```

To aid users in interpreting state, input, and parameter vector values, new methods for the [`separate_states`](@ref), [`separate_inputs`](@ref), and [`separate_parameters`](@ref) functions should be provided.  For the [`TypicalSection`](@ref) model, these new methods may be defined using the following block of code.

```julia
function separate_states(model::TypicalSection, x)

    return (h = x[1], theta = x[2], hdot = x[3], thetadot = x[4])
end

function separate_inputs(model::TypicalSection, y)

    return (L = y[1], M = y[2])
end

function separate_parameters(model::TypicalSection, p)

    return (kh = p[1], ktheta = p[1], m = p[1], Stheta = p[1], Itheta = p[1])
end
```

## Creating a New Model Coupling

In this section, we describe in detail how to construct a new model by coupling multiple existing models together.  To demonstrate this process, we also show how one might implement the [`Wagner`](@ref) and [`TypicalSection`](@ref) model coupling.

### Coupled Model Theory

To construct a coupled model, we first concatenate the governing differential equations, state variables, inputs, and parameters of multiple independent models
into a single system of equations.  This system of equations may be expressed as
```math
M(x,y,p,t)\dot{x} = f(x,y,p,t)
```
where
```math
x = \begin{bmatrix} x_1^T & x_2^T & \dots & x_N^T \end{bmatrix}^T \quad
y = \begin{bmatrix} y_1^T & y_2^T & \dots & y_N^T \end{bmatrix}^T \quad
p = \begin{bmatrix} p_1^T & p_2^T & \dots & p_N^T \end{bmatrix}^T \\
M(x, y, p, t) = \begin{bmatrix}
    M_1(x_1, y_1, p_1, t) & 0 & & 0 \\
    0 & M_2(x_2, y_2, p_2, t) & & 0 \\
    & & \ddots & \\
    0 & 0 & & M_N(x_N, y_N, p_N, t)
\end{bmatrix} \\
f(x, y, p, t) = \begin{bmatrix}
    f_1(x_1, y_1, p_1, t) \\
    f_2(x_2, y_2, p_2, t) \\
    \vdots \\
    f_N(x_N, y_N, p_N, t) \\
\end{bmatrix}
```

We then couple the models in this system of equations by allowing the inputs corresponding to each model to be defined as a function of the variables associated with the other models.  Specifically, we assume that the inputs corresponding to each model may be defined as
```math
y = g(x, p, t) - M_y(x, p, t) \dot{x}
```
where ``g(x, p, t)`` is a function which returns a vector and ``M_y(x, p, t)`` is a function which returns a matrix.  If we also assume that the state rates of the coupled model returned by ``f(x, y, p, t)`` are linearly dependent on the coupled model's inputs, at least for the non-zero rows of ``M_y(x, p, t)``, then the governing equations for the coupled system of equations may be expressed as
```math
\bar{M}(x, p, t) \dot{x} = \bar{f}(x, p, t)
```
where
```math
\bar{M}(x, p, t) = M(x, y, p, t) + \frac{\partial f}{\partial y} M_y(x, p, t)
```
and
```math
\bar{f}(x, p, t) = f(x, g(x, p, t), p, t)
```

For the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model, the coupled model's state variables, inputs, and parameters are
```math
x = \begin{bmatrix} \lambda_1 & \lambda_2 & h & \theta & \dot{h} & \dot{\theta} \end{bmatrix}^T \quad
y = \begin{bmatrix} u & v & \omega & L & M \end{bmatrix}^T \\
p = \begin{bmatrix} a & b & a_0 & \alpha_0 & k_h & k_\theta & m & S_\theta & I_\theta \end{bmatrix}^T \\
```
The [`Wagner`](@ref) model inputs may be defined as a function of the state variables of the typical section model.
```math
u = U_\infty \\
v = \dot{h} \\
\omega = \dot{\theta}
```
The [`TypicalSection`](@ref) model inputs may be defined using the lift and moment expressions corresponding to the [`Wagner`](@ref) model.
```math
\mathcal{L} = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u\omega - a b \dot{\omega} \right) \\
\mathcal{M} = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\omega + b \left( \frac{1}{8} - \frac{a}{2} \right) \dot{\omega} \right] + b \left(\frac{1}{2} + a \right) \mathcal{L}
```
where
```math
\alpha_\text{eff} = \left(\theta - \frac{v}{u} + \frac{b}{u} \left( \frac{1}{2} - a \right) \omega - \alpha_0 \right) \phi(0) + \frac{\lambda_1}{u} + \frac{\lambda_2}{u}
```

In order to use these relationships, two additional parameters must be introduced: the freestream velocity ``U_\infty`` and the freestream density ``\rho_\infty``.  These parameters are simply appended to the end of the parameter vector.

### Defining the Coupled Model

New coupled model's are defined by defining a new method for the [`couple_models`](@ref) function.  This function returns a tuple of models in the order in which their states, inputs, and parameters are concatenated.  A docstring should also be provided along with the method definition which defines the identities of any additional parameters, if used.

For the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model, a new coupled model may be defined using the following code block.

```julia
"""
    couple_models(aero::Wagner, stru::TypicalSection)

Create an aerostructural model using an unsteady aerodynamic model based on
Wagner's function and a two-degree of freedom typical section model.  This model
introduces the freestream velocity ``U_\\infty`` and air density ``\\rho_\\infty``
as additional parameters.
"""
couple_models(aero::Wagner, stru::TypicalSection) = (aero, stru)
```

### Defining the Coupled Model's Traits

The next step in defining a coupled model is to define the model's properties.  At a minimum, this requires defining new methods for the [`number_of_additional_parameters`](@ref) and [`coupling_inplaceness`](@ref AerostructuralDynamics.inplaceness) functions, though additional method definitions may be necessary.

The number of additional parameters introduced by the coupled model is specified by defining a new method for the [`number_of_additional_parameters`](@ref) function.  For out-of-place models, these methods must operate on the model types so that the length of the parameter vector is completely inferrable. For in-place models, this restriction is loosened and this methods may operate on model instances instead.

Whether the coupled model uses an in-place or out-of-place coupling function is defined by the [`coupling_inplaceness`](@ref AerostructuralDynamics.coupling_inplaceness) function, which operates on the model type. For performance reasons, in-place functions are generally preferred. The one exception is for coupled models with small numbers of state variables, in which case the preferred approach is to use static arrays with out-of-place functions.

The properties of the coupling function's mass matrix and/or state jacobian are defined by defining new methods for the [`coupling_mass_matrix_type`](@ref AerostructuralDynamics.mass_matrix_type) and/or [`coupling_state_jacobian_type`](@ref AerostructuralDynamics.state_jacobian_type) functions, respectively. By default, these properties assume their loosest possible definitions.

The properties of the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model may be defined using the following block of code.

```julia
number_of_additional_parameters(::Type{<:Wagner}, ::Type{TypicalSection}) = 2
coupling_inplaceness(::Type{<:Wagner}, ::Type{TypicalSection}) = OutOfPlace()
coupling_mass_matrix_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Linear()
coupling_state_jacobian_type(::Type{<:Wagner}, ::Type{TypicalSection}) = Nonlinear()
```

### Defining Coupling Function Methods

Once the properties of a coupled model have been defined, methods must be provided which define the values of the model's inputs. The portion of the inputs which is independent of the state rates is calculated using the [`get_coupling_inputs`](@ref) function for out-of-place coupling functions or the [`get_coupling_inputs!`](@ref) function for in-place coupling functions. For models with inputs that are also linearly dependent on the state rates, a new method must also be defined for the [`get_coupling_mass_matrix`](@ref AerostructuralDynamics.get_coupling_mass_matrix) function (or [`get_coupling_mass_matrix!`](@ref AerostructuralDynamics.get_coupling_mass_matrix!) function if the inputs are defined in-place). For constant mass matrices (coupling_mass_matrix_type(typeof.(models)...) == Constant()), this function should be defined without the x, p, and t arguments.

As an example, the expressions defining the inputs for the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model may be defined using the following block of code.

```julia
function get_coupling_inputs(aero::Wagner, stru::TypicalSection, s, p, t)
    # extract state variables
    λ1, λ2, h, θ, hdot, θdot = s
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    u = U
    v = U*θ + hdot
    ω = θdot
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L = tmp1*((v + d*ω - u*α0)*ϕ0 + λ1 + λ2) + tmp2*u/b*ω
    # moment at reference point
    M = -tmp2*u*ω + (b/2 + a*b)*L
    # return portion of inputs that is not dependent on the state rates
    return SVector(u, v, ω, L, M)
end

function get_coupling_mass_matrix(aero::Wagner, stru::TypicalSection, s, p, t)
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L_hddot = tmp/b
    L_θddot = -a*tmp
    # moment at reference point
    M_hddot = -tmp/2 + (b/2 + a*b)*L_hddot
    M_θddot = -tmp*(b/8 - a*b/2) + (b/2 + a*b)*L_θddot
    # construct submatrices
    Mda = @SMatrix [0 0; 0 0; 0 0]
    Mds = @SMatrix [0 0 0 0; 0 0 0 0; 0 0 0 0]
    Mra = @SMatrix [0 0; 0 0]
    Mrs = @SMatrix [0 0 -L_hddot -L_θddot; 0 0 -M_hddot -M_θddot]
    # assemble mass matrix
    return [Mda Mds; Mra Mrs]
end
```

### Defining Methods for the Coupling Function's Jacobians

Unless otherwise specified, the jacobian of the coupling function of a given coupled model with respect to the state variables is calculated when necessary using forward automatic differentiation (enabled by the ForwardDiff package). While this approach for computing the jacobians is convenient and exact, alternative methods for computing jacobians may be more computationally efficient. To manually define the jacobian of the coupling function with respect to the state variables, a new method for the [`get_coupling_state_jacobian`](@ref AerostructuralDynamics.get_coupling_state_jacobian) (or [`get_coupling_state_jacobian!`](@ref AerostructuralDynamics.get_coupling_state_jacobian) for inputs which are defined in-place) may be defined.

For the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model, this jacobian may be defined analytically using the following block of code.

```julia
function get_coupling_state_jacobian(aero::Wagner, stru::TypicalSection, u, p, t) where {N,TF,SV,SA}
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # extract model constants
    C1 = aero.C1
    C2 = aero.C2
    # local freestream velocity components
    v_θ = U
    v_hdot = 1
    ω_θdot = 1
    # calculate loads
    r_λ = wagner_loads_λ(a, b, ρ, a0, U)
    L_h, M_h = wagner_loads_h()
    L_θ, M_θ = wagner_loads_θ(a, b, ρ, a0, C1, C2, U)
    L_hdot, M_hdot = wagner_loads_v(a, b, ρ, a0, C1, C2, U)
    L_θdot, M_θdot = wagner_loads_ω(a, b, ρ, a0, C1, C2, U)
    # compute jacobian sub-matrices
    Jda = @SMatrix [0 0; 0 0; 0 0]
    Jds = @SMatrix [0 0 0 0; 0 v_θ v_hdot 0; 0 0 0 ω_θdot]
    Jra = r_λ
    Jrs = @SMatrix [L_h L_θ L_hdot L_θdot; M_h M_θ M_hdot M_θdot]
    # return jacobian
    return [Jda Jds; Jra Jrs]
end

function wagner_loads_λ(a, b, ρ, a0, u)
    tmp1 = a0*ρ*u*b
    tmp2 = (b/2 + a*b)*tmp1
    return @SMatrix [tmp1 tmp1; tmp2 tmp2]
end

wagner_loads_h() = SVector(0, 0)

function wagner_loads_θ(a, b, ρ, a0, C1, C2, u)
    ϕ0 = 1 - C1 - C2
    L_θ = a0*ρ*b*u^2*ϕ0
    M_θ = (b/2 + a*b)*L_θ
    return SVector(L_θ, M_θ)
end

function wagner_loads_v(a, b, ρ, a0, C1, C2, u)
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_v = a0*ρ*u*b*ϕ0
    # moment at reference point
    M_v = (b/2 + a*b)*L_v

    return SVector(L_v, M_v)
end

function wagner_loads_ω(a, b, ρ, a0, C1, C2, u)
    # circulatory load factor
    tmp1 = a0*ρ*u*b
    # non-circulatory load factor
    tmp2 = pi*ρ*b^3
    # constant based on geometry
    d = b/2 - a*b
    # Wagner's function at t = 0.0
    ϕ0 = 1 - C1 - C2
    # lift at reference point
    L_ω = tmp1*d*ϕ0 + tmp2*u/b
    # moment at reference point
    M_ω = -tmp2*u + (b/2 + a*b)*L_ω

    return SVector(L_ω, M_ω)
end
```

### Defining Methods for Unit Testing

In order to test whether the provided mass matrices are correct for a given coupled model, a new method for [`get_coupling_inputs_using_state_rates`](@ref AerostructuralDynamics.get_inputs_using_state_rates) (which defines the portion of the inputs that are dependent on the state rates) must be provided.  Since this function is used for testing, there is no in-place version of this function.

For the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model, this function could be defined as follows

```julia
function get_inputs_using_state_rates(aero::Wagner, stru::TypicalSection,
    ds, s, p, t)
    # extract state rates
    dλ1, dλ2, dh, dθ, dhdot, dθdot = ds
    # extract parameters
    a, b, a0, α0, kh, kθ, m, Sθ, Iθ, U, ρ = p
    # local freestream velocity components
    udot = 0
    vdot = dhdot
    ωdot = dθdot
    # non-circulatory load factor
    tmp = pi*ρ*b^3
    # lift at reference point
    L = tmp*(vdot/b - a*ωdot)
    # moment at reference point
    M = -tmp*(vdot/2 + b*(1/8 - a/2)*ωdot) + (b/2 + a*b)*L
    # return inputs
    return SVector(0, 0, 0, L, M)
end
```

### Defining Convenience Methods

To aid users in defining state, input, and parameter vectors, new methods for the [`set_states!`](@ref), [`set_inputs!`](@ref), and [`set_parameters!`](@ref) functions should be provided.  For the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model, these new methods may be defined as follows

```julia
function set_states!(x, aero::Wagner, stru::TypicalSection; lambda = [0.0, 0.0],
    h = 0.0, theta = 0.0, hdot = 0.0, thetadot = 0.0)

    models = (aero, stru)
    xs = view.(x, state_indices.(models))

    set_states!(xs[1], aero; lambda)
    set_states!(xs[2], stru; h, theta, hdot, thetadot)

    return x
end

function set_inputs!(y, aero::Wagner, stru::TypicalSection; u=1.0, v=0.0,
    omega=0.0, L=0.0, M=0.0)

    models = (aero, stru)
    ys = view.(y, input_indices.(models))

    set_inputs!(ys[1], aero; u, v, omega)
    set_inputs!(ys[2], stru; L, M)

    return y
end

function set_parameters!(p, aero::Wagner, stru::TypicalSection; a=0.0, b=0.5,
    a0=2*pi, alpha0=0.0, kh, ktheta, m, Stheta, Itheta)

    models = (aero, stru)
    ps = view.(p, input_indices.(models))

    set_inputs!(ps[1], aero; a, b, a0, alpha0)
    set_inputs!(ps[2], stru; kh, ktheta, m, Stheta, Itheta)

    return y
end
```

To aid users in interpreting state, input, and parameter vector values, new methods for the [`separate_states`](@ref), [`separate_inputs`](@ref), and [`separate_parameters`](@ref) functions should be provided.  For the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model, these new methods may be defined using the following block of code.

```julia
function separate_states(model::TypicalSection, x)

    return (h = x[1], theta = x[2], hdot = x[3], thetadot = x[4])
end

function separate_inputs(model::TypicalSection, y)

    return (L = y[1], M = y[2])
end

function separate_parameters(model::TypicalSection, p)

    return (kh = p[1], ktheta = p[1], m = p[1], Stheta = p[1], Itheta = p[1])
end
```

### Defining Visualization Methods

To aid users in visualizing solution geometry new plot recipes

defining state, input, and parameter vectors, new methods for the [`set_states!`](@ref), [`set_inputs!`](@ref), and [`set_parameters!`](@ref) functions should be provided.  For the [`Wagner`](@ref) model coupled with the [`TypicalSection`](@ref) model, these new methods are defined as follows

### Model Ordering

In general, we suggest that the following model order is used when constructing coupled models.

1. Aerodynamics Model(s)
2. Structural Model(s)
3. Rigid Body Model(s)
4. Control Surface Model(s)
5. Controller Model(s)

## Avoiding Mass Matrices

In order to take advantage of as many features of the DifferentialEquations package as possible (including local sensitivity analysis), at this point in time we recommend that the governing differential equations for models be reformulated to avoid using mass matrices whenever possible.
