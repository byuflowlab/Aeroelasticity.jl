# Developer's Guide

In this guide we demonstrate by example how to define aerodynamic and structural models and couple them together.

```@contents
Pages = ["developer.md"]
Depth = 3
```

```@setup developer
using AerostructuralDynamics, StaticArrays, LinearAlgebra
```

```@setup developer-combined
using AerostructuralDynamics, StaticArrays, LinearAlgebra #hide
```

## Defining a Structural Model

For demonstrating how to create a new structural model, we use a typical section model with two degrees of freedom, as shown in the following figure.

![](typical-section.svg)

The equations of motion for this model are:

```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M_{\frac{1}{4}} + b \left( \frac{1}{2} + a \right) L
```

where ``a`` is the normalized distance from the semichord to the reference point, ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_P`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M_\frac{1}{4}`` is the quarter-chord moment per unit span.

### Theory

Structural state variables in this package satisfy the ordinary differential equation (or differential algebraic equation in mass matrix form)

```math
M_s(q,r,p_s,t)\dot{q} = f_s(q,r,p_s,t)
```

where ``M_s`` is a function which defines the structural mass matrix, ``f_s`` is a function which defines the mass matrix multiplied structural state rates, ``q`` is a vector of structural states, ``r`` is a vector of aerodynamic loads, ``p_s`` is a vector of structural parameters, and ``t`` is the current time.

The equations of motion for the typical section model when expressed in the form expected by this package are

```math
M \dot{q} = K q + D r
```

where

```math
q = \begin{bmatrix} h & \theta & \dot{h} & \dot{\theta} \end{bmatrix}^T \quad
r = \begin{bmatrix} L & M_\frac{1}{4} \end{bmatrix}^T
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
b \left( \frac{1}{2} + a \right) & 1
\end{bmatrix}
```

### Defining a New Type

We start creating our model by defining a new type.  

```@example developer
"""
    TypicalSection <: StructuralModel

Typical section structural model with state variables ``q = \\begin{bmatrix} h &
θ & \\dot{h} & \\dot{\\theta} \\end{bmatrix}^T``, structural parameters ``p_s =
\\begin{bmatrix} a & b & k_h & k_\\theta & m & x_\\theta & I_P \\end{bmatrix}^T``,
and aerodynamic loads ``r = \\begin{bmatrix} L & M_\\frac{1}{4} \\end{bmatrix}^T``
"""
struct TypicalSection <: StructuralModel end

nothing #hide
```

Model parameters may be stored either as fields of the newly defined type or passed directly to the solver.  In this case, we choose to define all of the structural parameters as elements of the parameter vector.  In general, model constants should be stored as fields of the struct, whereas parameters that may change should be passed to the solver through the parameter vector.

Note that the state, parameter, and aerodynamic load variable identities should be documented in the docstring of the type definition since the type definition provides the primary (and sometimes only) source of documentation for a given model.

### Defining Model Properties

We now need to define a few model properties.

To indicate the number of state variables in our model, we define a new method for the [`number_of_states`](@ref) function.  In our case, we have four structural state variables ``q = \begin{bmatrix} h & \theta & \dot{h} & \dot{\theta} \end{bmatrix}^T ``.  This method must be defined for all models.

```@example developer
number_of_states(::TypicalSection) = 4
nothing #hide
```

To indicate the number of aerodynamic load inputs, we define a new method for the [`number_of_inputs`](@ref) function.  In our case, we have two aerodynamic loads ``r = \begin{bmatrix} L & M_\frac{1}{4} \end{bmatrix}^T``.  This method must be defined for all models.

```@example developer
number_of_inputs(::TypicalSection) = 2
nothing #hide
```

To indicate whether we plan to use in-place or out-of-place functions, we define a new method for the [`isinplace`](@ref) function.  In general, for performance reasons, in-place functions are preferred.  The one exception is for models with small amounts of state variables, in which case the preferred approach is to use static arrays with out-of-place functions.  For this model, we use the latter approach.  This method must be defined for all models.

```@example developer
isinplace(::TypicalSection) = false
nothing #hide
```

To indicate whether we plan to use a mass matrix, we define a new method for the [`has_mass_matrix`](@ref) function.  This method must be defined for all models.

```@example developer
has_mass_matrix(::TypicalSection) = true
nothing #hide
```

To indicate whether the mass matrix is constant (not dependent on the state variables, aerodynamic loads, parameters, or time), we define a new method for the [`constant_mass_matrix`](@ref) function.  For the typical section model, the mass matrix is not constant since it is dependent on the model parameters.  This function definition is required if `has_mass_matrix(model) == true`

```@example developer
constant_mass_matrix(::TypicalSection) = false
nothing #hide
```

To indicate whether the state rates are linearly dependent on the aerodynamic loads, or in other words whether the governing structural equations can be expressed as  
```math
M_s(q,p_s,t)\dot{q} = f_s(q,p_s,t) + D(q,p_s,t)r
```
we define a new method for the [`linear_input_dependence`](@ref) function.  Having a linear load dependence allows for greater flexibility when coupling the structural model with an aerodynamic model.  The default return value from this function is `false`.

```@example developer
linear_input_dependence(::TypicalSection) = true
nothing #hide
```

To indicate whether the jacobian of the right hand side of the rate equations with respect to the state variables is manually defined, we define a new method for the
[`defined_state_jacobian`](@ref) function.  The default return value from this function is `false`.

```@example developer
defined_state_jacobian(::TypicalSection) = true
nothing #hide
```

To indicate whether the jacobian of the right hand side of the rate equations with respect to the aerodynamic loads is manually defined, we define a new method for the
[`defined_input_jacobian`](@ref) function.  The default return value from this function is `false`.

```@example developer
defined_input_jacobian(::TypicalSection) = true
nothing #hide
```

To indicate whether the jacobian of the right hand side of the rate equations with respect to the aerodynamic loads is constant (not dependent on the state variables, aerodynamic loads, parameters, or time), we define a new method for the
[`constant_input_jacobian`](@ref) function.  The default return value from this function is `false`.

```@example developer
constant_input_jacobian(::TypicalSection) = false
nothing #hide
```

### Mass Matrix Equation

For out-of-place models, the mass matrix is calculated using the [`get_mass_matrix`](@ref) function.  For in-place models the mass matrix is calculated using the [`get_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `q`, `r`, `p`, and `t` arguments.  

```@example developer
function get_mass_matrix(::TypicalSection, q, r, p, t)    
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # update mass matrix
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
end
nothing #hide
```

### State Rate Equation

The right hand side of the governing structural differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models and [`get_rates!`](@ref) function for in-place models.  

```@example developer
function get_rates(::TypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q    
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # calculate state rates
    return SVector(hdot, θdot, -kh*h - L, -kθ*θ + M + b*(1/2+a)*L)
end

nothing #hide
```

### State Rate Jacobian

The jacobian of the right hand side of the governing structural equations with respect to the state variables is calculated using the [`get_state_jacobian`](@ref) function for out-of-place models and [`get_state_jacobian!`](@ref) function for in-place models.

```@example developer
function get_state_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

nothing #hide
```

Note that this function is only used if `defined_state_jacobian(model) == true`.

### Aerodynamic Load Jacobian

The jacobian of the right hand side of the governing structural equations with respect to the aerodynamic loads is defined using the [`get_load_jacobian`](@ref) function.  If the jacobian of the right hand side of the governing structural equations with respect to the aerodynamic loads is constant `constant_input_jacobian(stru)==true`, this function is called without the `q`, `r`, `p`, and `t` arguments.  There is no out-of-place form for this function, however, it may be constructed as a either a linear map (if large) or static array (if small) in order to avoid allocations.  

```@example developer
function get_load_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0; 0 0; -1 0; b*(1/2 + a) 1]
end

nothing #hide
```

Note that this function is only used if `defined_input_jacobian(model) == true`.

### Structural Model Code

Putting it all together, a complete representation of our typical section model for use with this package may be defined using the following block of code.

```@example developer-combined
"""
    TypicalSection <: StructuralModel

Typical section structural model with state variables ``q = \\begin{bmatrix} h &
θ & \\dot{h} & \\dot{\\theta} \\end{bmatrix}^T``, structural parameters ``p_s =
\\begin{bmatrix} a & b & k_h & k_\\theta & m & x_\\theta & I_P \\end{bmatrix}^T``,
and aerodynamic loads ``r = \\begin{bmatrix} L & M_\\frac{1}{4} \\end{bmatrix}^T``
"""
struct TypicalSection <: StructuralModel end

number_of_states(::TypicalSection) = 4
number_of_inputs(::TypicalSection) = 2
isinplace(::TypicalSection) = false
has_mass_matrix(::TypicalSection) = true
constant_mass_matrix(::TypicalSection) = false
linear_input_dependence(::TypicalSection) = true
defined_state_jacobian(::TypicalSection) = true
defined_input_jacobian(::TypicalSection) = true
constant_input_jacobian(::TypicalSection) = false

function get_mass_matrix(::TypicalSection, q, r, p, t)    
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # update mass matrix
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
end

function get_rates(::TypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q    
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # calculate state rates
    return SVector(hdot, θdot, -kh*h - L, -kθ*θ + M + b*(1/2+a)*L)
end

function get_state_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

function get_load_jacobian(::TypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0; 0 0; -1 0; b*(1/2 + a) 1]
end

nothing #hide
```

## Defining an Aerodynamic Model

For demonstrating how to create a new aerodynamic model, we use Peter's finite state aerodynamic model.

The governing differential equation for the aerodynamic states is

```math
\bar{A}\dot{\lambda} + \frac{U}{b}\lambda = \bar{c}\left[ \ddot{h} + U\dot{\theta} + b \left(\frac{1}{2} - a\right) \ddot{\theta} \right]
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

When aerodynamic models are used in isolation, the structural deflections are appended to the end of the parameter vector.  When coupled with a structural model, the structural deflections are calculated in the format expected by the aerodynamic model as a function of the structural and aerodynamic states.  The code for defining an aerodynamic model is independent of the choice of structural model, so the same set of functions is used when an aerodynamic model is considered in isolation and when it is coupled with a structural model  

### Theory

Aerodynamic state variables in this package satisfy the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M_{a}(λ,d,p_a,t)\dot{\lambda} = f_a(λ,d,p_a,t)
```

where ``M_a`` is a function which defines the aerodynamic mass matrix, ``f_a`` is a function which defines the mass matrix multiplied aerodynamic state rates, ``\lambda`` is a vector of structural states, ``d`` is a vector of structural deflections, ``p_a`` is a vector of aerodynamic parameters, and ``t`` is the current time.

The equations of motion for Peter's finite state model, when expressed in the form expected by this package are:

```math
\bar{A}\dot{\lambda} = -\frac{U}{b}\lambda +
E d
```
where
```math
E = \bar{c} \begin{bmatrix}
U \\
1 \\
b \left(\frac{1}{2} - a\right)
\end{bmatrix}
\quad
d = \begin{bmatrix}
\dot{\theta} \\
\ddot{h} \\
\ddot{\theta}
\end{bmatrix}
```

### Defining a New Type

We start creating our aerodynamic model by defining a new type.  

```@example developer
"""
    PetersFiniteState{N, TF} <: AerodynamicModel

Peter's finite state model with `N` aerodynamic states, aerodynamic parameters
``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``, and structural
deflections ``d = \\begin{bmatrix} \\dot{\\theta} & \\ddot{h} &
\\ddot{\\theta}\\end{bmatrix}^T``
"""
struct PetersFiniteState{N, TF} <: StructuralModel
    A::SMatrix{N,N,TF}
    b::SVector{N,TF}
    c::SVector{N,TF}
end

nothing #hide
```

Here `N` is the number of aerodynamic state variables and `TF` is the floating point type used to represent the constant matrices/vectors ``\bar{A}``, ``\bar{b}``, and ``\bar{c}``.

For convenience, we create a constructor which initializes matrix ``\bar{A}`` and vectors ``\bar{b}`` and ``\bar{c}`` given the number of aerodynamic state variables and floating point type.

```@example developer
"""
    PetersFiniteState{N,TF}()

Initialize an object of type `PetersFiniteState` which has `N` aerodynamic
degrees of freedom.
"""
function PetersFiniteState{N,TF}() where {N,TF}

    b = zeros(TF, N)
    for i = 1:N
        b[i] = (-1)^(n-1)*factorial(N + n)/factorial(N - n)*
            1/factorial(n)^2
    end
    b[N] += (-1)^N

    c = zeros(TF, N)
    c .= 2/n

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[m, n] = n/2
    end
    for m in 2:N
        n = m - 1
        D[m, n] = -n/2
    end

    A = D + d*b' + c*d' + 1//2*c*b'

    return PetersFiniteState{N,TF}(A,b,c)
end

nothing #hide
```

### Defining Model Properties

We now need to define a few model properties.

To indicate the number of state variables in our model, we define a new method for the [`number_of_states`](@ref) function.  In our case, we have an arbitrary number of aerodynamic states ``\lambda``, though typically 3-10 aerodynamic states are used with Peter's finite state model.  This method must be defined for all models.

```@example developer
number_of_states(::PetersFiniteState{N,TF}) where {N,TF} = N
nothing #hide
```

To indicate the number of structural deflection inputs, we define a new method for the [`number_of_inputs`](@ref) function.  In our case, we have three structural deflection variables ``d = \begin{bmatrix} \dot{\theta} & \ddot{h} & \ddot{\theta} \end{bmatrix}^T``.  This method must be defined for all models.

```@example developer
number_of_inputs(::PetersFiniteState) = 3
nothing #hide
```

To indicate whether we plan to use in-place or out-of-place functions, we define a new method for the [`isinplace`](@ref) function.  This method must be defined for all models.

```@example developer
isinplace(::PetersFiniteState) = false
nothing #hide
```

To indicate whether we plan to use a mass matrix, we define a new method for the [`has_mass_matrix`](@ref) function.  This method must be defined for all models.

```@example developer
has_mass_matrix(::PetersFiniteState) = true
nothing #hide
```

To indicate whether the mass matrix is constant (not dependent on the state variables, structural deflections, parameters, or time), we define a new method for the [`constant_mass_matrix`](@ref) function.  For Peter's finite state model, the mass matrix is constant since it depends only on the number of aerodynamic states.  This function definition is required if `has_mass_matrix(model) == true`

```@example developer
constant_mass_matrix(::PetersFiniteState) = true
nothing #hide
```

To indicate whether the state rates are linearly dependent on the structural deflections, or in other words, whether the governing aerodynamic equations may be expressed as
```math
M_a(\lambda,p_a,t)\dot{\lambda} = f_a(\lambda,p_a,t) + E(\lambda,p_a,t)d
```
we define a new method for the [`linear_input_dependence`](@ref) function.  Having a linear deflection dependence allows for greater flexibility when coupling the aerodynamic model with a structural model.  The default return value from this function is `false`.

```@example developer
linear_input_dependence(::PetersFiniteState) = true
nothing #hide
```

To indicate whether the jacobian of the right hand side of the rate equations with respect to the state variables is manually defined, we define a new method for the
[`defined_state_jacobian`](@ref) function.  The default return value from this function is `false`.

```@example developer
defined_state_jacobian(::PetersFiniteState) = true
nothing #hide
```

To indicate whether the jacobian of the right hand side of the rate equations with respect to the structural deflections is manually defined, we define a new method for the
[`defined_input_jacobian`](@ref) function.  The default return value from this function is `false`.

```@example developer
defined_input_jacobian(::PetersFiniteState) = true
nothing #hide
```

To indicate whether the jacobian of the right hand side of the rate equations with respect to the structural deflections is constant (not dependent on the state variables, aerodynamic loads, parameters, or time), we define a new method for the
[`constant_input_jacobian`](@ref) function.  The default return value from this function is `false`.

```@example developer
constant_input_jacobian(::PetersFiniteState) = false
nothing #hide
```

### Mass Matrix Equation

For out-of-place models, the mass matrix is calculated using the [`get_mass_matrix`](@ref) function.  For in-place models the mass matrix is calculated using the [`get_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `λ`, `d`, `p`, and `t` arguments.  

```@example developer
get_mass_matrix(model::PetersFiniteState) =  model.A
nothing #hide
```

### State Rate Equation

The right hand side of the governing aerodynamic differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models and [`get_rates!`](@ref) function for in-place models.  

```@example developer
function get_rates(model::PetersFiniteState, λ, d, p, t)
    # extract structural deflections
    θdot, hddot, θddot = d
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return cbar*(hddot + U*θdot + b*(1/2-a)*θddot) - U/b*λ
end

nothing #hide
```

### State Rate Jacobian

The jacobian of the right hand side of the governing aerodynamic equations with respect to the state variables is calculated using the [`get_state_jacobian`](@ref) function for out-of-place models and [`get_state_jacobian!`](@ref) function for in-place models.

```@example developer
function get_state_jacobian(model::PetersFiniteState, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # jacobian with respect to aerodynamic states
    return -U/b*Diagonal(one.(model.c))
end

nothing #hide
```

Note that this function is only used if `defined_state_jacobian(model) == true`.

### Structural Deflection Jacobian

The jacobian of the right hand side of the governing aerodynamic equations with respect to the structural deflections is defined using the [`get_deflection_jacobian`](@ref) function.  If the jacobian of the right hand side of the governing aerodynamic equations with respect to the structural deflections is constant (`constant_input_jacobian(aero)==true`), this function is called without the `λ`, `d`, `p`, and `t` arguments. There is no out-of-place form for this function, however, it may be constructed as a either a linear map (if large) or static array (if small) in order to avoid allocations.

```@example developer
function get_deflection_jacobian(::PetersFiniteState, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return hcat(U*cbar, cbar, b*(1/2-a)*cbar)
end

nothing #hide
```

Note that this function is only used if `defined_input_jacobian(model) == true`.

### Aerodynamic Model Code

Putting it all together, a complete representation Peter's finite state aerodynamic model may be defined using the following block of code.

```@example developer-combined

"""
    PetersFiniteState{N, TF} <: AerodynamicModel

Peter's finite state model with `N` aerodynamic states, aerodynamic parameters
``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``, and structural
deflections ``d = \\begin{bmatrix} \\dot{\\theta} & \\ddot{h} &
\\ddot{\\theta}\\end{bmatrix}^T``
"""
struct PetersFiniteState{N, TF} <: StructuralModel
    A::SMatrix{N,N,TF}
    b::SVector{N,TF}
    c::SVector{N,TF}
end

"""
    PetersFiniteState{N,TF}()

Initialize an object of type `PetersFiniteState` which has `N` aerodynamic
degrees of freedom.
"""
function PetersFiniteState{N,TF}() where {N,TF}

    b = zeros(TF, N)
    for i = 1:N
        b[i] = (-1)^(n-1)*factorial(N + n)/factorial(N - n)*
            1/factorial(n)^2
    end
    b[N] += (-1)^N

    c = zeros(TF, N)
    c .= 2/n

    d = zeros(TF, N)
    d[1] = 1/2

    D = zeros(TF, N, N)
    for m in 1:N-1
        n = m + 1
        D[m, n] = n/2
    end
    for m in 2:N
        n = m - 1
        D[m, n] = -n/2
    end

    A = D + d*b' + c*d' + 1//2*c*b'

    return PetersFiniteState{N,TF}(A,b,c)
end

number_of_states(::PetersFiniteState{N,TF}) where {N,TF} = N
number_of_inputs(::PetersFiniteState) = 3
isinplace(::PetersFiniteState) = false
has_mass_matrix(::PetersFiniteState) = true
constant_mass_matrix(::PetersFiniteState) = true
linear_input_dependence(::PetersFiniteState) = true
defined_state_jacobian(::PetersFiniteState) = true
defined_input_jacobian(::PetersFiniteState) = true
constant_input_jacobian(::PetersFiniteState) = false

get_mass_matrix(model::PetersFiniteState) =  model.A

function get_rates(model::PetersFiniteState, λ, d, p, t)
    # extract structural deflections
    θdot, hddot, θddot = d
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return cbar*(hddot + U*θdot + b*(1/2-a)*θddot) - U/b*λ
end

function get_state_jacobian(model::PetersFiniteState, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # jacobian with respect to aerodynamic states
    return -U/b*Diagonal(one.(model.c))
end

function get_deflection_jacobian(::PetersFiniteState, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return hcat(U*cbar, cbar, b*(1/2-a)*cbar)
end

nothing #hide
```

## Coupling Aerodynamic and Structural Models

For demonstrating how to couple aerodynamic and structural models together, we use the typical section structural model and Peter's finite state aerodynamic model.  To couple these two models together, we need to define the model inputs (aerodynamic loads and structural deflections) as functions of the states, parameters, and time.

### Theory

We assume the aerodynamic loads and structural deflections may be expressed as a function of the aerodynamic and structural state variables and parameters, as well as the current time.

```math
r = f_r(u,p,t) \quad d = f_d(u,p,t) \\
```
where
```math
u = \begin{bmatrix} q & \lambda \end{bmatrix}^T \quad
p = \begin{bmatrix} p_s & p_a \end{bmatrix}^T \quad
```

If the state rates of a given aerodynamic model are linearly dependent on the structural deflections (`linear_input_dependence(aero) == true`), we can expand the expression which defines the structural deflections to
```math
d = f_d(u,p,t) - M_{d s}(u,p,t) \dot{q} -
M_{d a}(u,p,t) \dot{\lambda}
```
where ``M_{d s}`` is a function which defines the (negative) jacobian of the structural deflections with respect to the structural state rates and ``M_{d a}`` is a function which defines the (negative) jacobian of the structural deflections with respect to the aerodynamic state rates.  ``f_d`` is a function which defines the portion of the structural deflections which is independent of the structural and aerodynamic state rates.

If the state rates of a given structural model are linearly dependent on the aerodynamic loads (`linear_input_dependence(stru) == true`), we can expand the expression which defines the aerodynamic loads to
```math
r = f_r(u,p,t) - M_{r s}(u,p,t) \dot{q} -
M_{r a}(u,p,t) \dot{\lambda}
```
where ``M_{r s}`` is a function which defines the (negative) jacobian of the aerodynamic loads with respect to the structural state rates and ``M_{r a}`` is a function which defines the (negative) jacobian of the aerodynamic loads with respect to the aerodynamic state rates.  ``f_r`` is a function which defines the portion of the aerodynamic loads which is independent of the structural and aerodynamic state rates.

In the most general case, the coupled system of equations may be defined as
```math
M_u(u,p,t) \dot{u} = f_u(u,p,t)
```
where
```math
f_u = \begin{bmatrix} f_s \\f_a \end{bmatrix} \quad
M_{u} = \begin{bmatrix} M_s & 0 \\ 0 & M_a \end{bmatrix} + \begin{bmatrix} \frac{\partial f_s}{\partial r} & 0 \\ 0 & \frac{\partial f_a}{\partial d} \end{bmatrix} \begin{bmatrix} M_{rs} & M_{ra} \\ M_{ds} & M_{da} \end{bmatrix} \quad
```
The associated jacobian is
```math
\frac{\partial f_u}{\partial u} = \begin{bmatrix} \frac{\partial f_s}{\partial q} & 0 \\ 0 & \frac{\partial f_a}{\partial \lambda} \end{bmatrix} + \begin{bmatrix} \frac{\partial f_s}{\partial r} & 0 \\ 0 & \frac{\partial f_a}{\partial d} \end{bmatrix} \begin{bmatrix} \frac{\partial f_r}{\partial q} & \frac{\partial f_r}{\partial \lambda} \\ \frac{\partial f_d}{\partial q} & \frac{\partial f_d}{\partial \lambda} \end{bmatrix}
```

If we introduce the combined input function
```math
y = f_y(u, p, t) - M_y(u, p, t) \dot{u}
```
where
```math
y = \begin{bmatrix} r \\ d \end{bmatrix} \quad f_y = \begin{bmatrix} f_r \\ f_d \end{bmatrix} \quad M_y = \begin{bmatrix} M_{rs} & M_{ra} \\ M_{ds} & M_{da} \end{bmatrix}
```
the mass matrix and jacobian expressions may be shortened to
```math
M_{u} = \begin{bmatrix} M_s & 0 \\ 0 & M_a \end{bmatrix} + \begin{bmatrix} \frac{\partial f_s}{\partial r} & 0 \\ 0 & \frac{\partial f_a}{\partial d} \end{bmatrix} M_y
```
```math
\frac{\partial f_u}{\partial u} =
\begin{bmatrix}
\frac{\partial f_s}{\partial q} & 0 \\
0 & \frac{\partial f_a}{\partial \lambda}
\end{bmatrix} +
\begin{bmatrix}
\frac{\partial f_s}{\partial r} & 0 \\
0 & \frac{\partial f_a}{\partial d}
\end{bmatrix} \frac{\partial f_y}{\partial u}
```

We have implemented all of the functions in these expressions in the individual model definitions with the exception of those associated with the combined input function, so all that is required to couple aerodynamic and structural models together is to define the functions associated with the combined input function.

#### Aerodynamic Loads

The aerodynamic loads expected by the typical section model are the lift and quarter-chord moment.  The lift and quarter-chord moment, as calculated using Peter's finite state model are
```math
L = \pi \rho_\infty b^2 \left( \ddot{h} + U \dot{\theta} - b a \ddot{\theta} \right) + 2 \pi \rho_\infty U b \left[ h + U \theta + b \left( \frac{1}{2} - a \right) \dot{\theta} - \lambda_0 \right]
```
```math
M_\frac{1}{4} = - \pi \rho_\infty b^3 \left[ \frac{1}{2} \ddot{h} + U \dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right]
```
where ``U`` is the freestream velocity, ``\rho_\infty`` is the freestream air density, and ``\lambda_0`` is the induced flow velocity.  The induced flow velocity may be approximated as a function of the aerodynamic states
```math
\lambda_0 \approx \frac{1}{2} \sum\limits_{n=1}^{N}\bar{b}_n \lambda_n
```

The aerodynamic loads calculated by Peter's finite state model when coupled with the typical section model may be expressed in the form expected by this package as
```math
r = J_{rs} q + J_{ra} \lambda
- M_{rs} \dot{q}
```

where

```math
r =\begin{bmatrix}
L &
M_\frac{1}{4}
\end{bmatrix}^T
\quad
q = \begin{bmatrix}
h &
\theta &
\dot{h} &
\dot{\theta}
\end{bmatrix}^T
\quad
\lambda = \begin{bmatrix}
\lambda_1 &
\lambda_2 &
... &
\lambda_N
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

#### Structural Deflections

The structural deflections expected by Peter's finite state model are ``d = \begin{bmatrix} \dot{\theta} & \ddot{h} & \ddot{\theta} \end{bmatrix}^T``.  These structural deflections correspond to a subset of the structural states and corresponding rates.

The structural deflections may be expressed in the form expected by this package as
```math
d = J_{ds} q - M_{ds} \dot{q}
```
where
```math
J_{ds} = \begin{bmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \quad
M_{ds} = \begin{bmatrix} 0 & 0 & 0 & 0 \\ 0 & 0 & -1 & 0 \\ 0 & 0 & 0 & -1 \end{bmatrix}
```

#### Combined Input Function

Upon combining the expressions for the aerodynamic loads and structural deflections, we obtain the following expression for the combined input function.

```math
y = J_{y} y - M_{y} \dot{u}
```
where
```math
J_{y} = \begin{bmatrix} J_{rs} & J_{ra} \\ J_{ds} & 0 \end{bmatrix} \quad
M_{y} = \begin{bmatrix} M_{rs} & 0 \\ M_{ds} & 0 \end{bmatrix}
```

### Defining Input Function Properties

Before we define the combined input function (and its associated functions) we need to define a few of its properties.

To define whether the combined input function should use an in-place or out-of-place format, we define a new method for the [`inplace_input`](@ref) function.

```@example developer
inplace_input(::PetersFiniteState, ::TypicalSection) = false
nothing #hide
```

To indicate whether the combined input function involves a mass matrix term, we define a new method for the [`has_input_mass_matrix`](@ref) function.

```@example developer
has_input_mass_matrix(::PetersFiniteState, ::TypicalSection) = true
nothing #hide
```

To indicate whether the combined input function mass matrix is constant, we define a new method for the [`constant_input_mass_matrix`](@ref) function.  This function is required if `has_input_mass_matrix(aero, stru) == true`

```@example developer
constant_input_mass_matrix(::PetersFiniteState, ::TypicalSection) = false
nothing #hide
```

To indicate whether the jacobian of the combined input function ``f_y`` with respect to the state variables is defined, we define a new method for the [`defined_input_state_jacobian`] function.  The default return value from this function is `false`.

```@example developer
defined_input_state_jacobian(::PetersFiniteState, ::TypicalSection) = true
nothing #hide
```

### Input Mass Matrix Equation

For out-of-place combined input functions, the mass matrix is calculated using the [`get_input_mass_matrix`](@ref) function.  For in-place combined input functions, the mass matrix is calculated using the [`get_input_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `u`, `p`, and `t` arguments.

```@example developer

function get_input_mass_matrix(stru::TypicalSection, aero::PetersFiniteState, u, p, t)
    # extract model constants
    cbar = model.c
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # create zero row vector with length nλ
    zλ = zero(aero.c')
    # construct submatrices
    Mrs = pi*ρ*b^2 * @SMatrix [0 0 -1 a*b; 0 0 b/2 b^2*(1/8-a/2)]
    Mra = @SMatrix vcat(zλ, zλ)
    Mds = @SMatrix [0 0 0 0; 0 0 -1 0; 0 0 0 -1]
    Mda = @SMatrix vcat(zλ, zλ, zλ, zλ)
    # assemble mass matrix
    return vcat(hcat(Mrs, Mra), hcat(Mds, Mda))
end

nothing #hide
```

### Input Equation

The portion of the inputs which is independent of the state rates is calculated using the [`get_inputs`](@ref) function for out-of-place aerodynamic load calculations and [`get_inputs!`](@ref) function for in-place aerodynamic load calculations.

```@example developer

function get_inputs(stru::TypicalSection, aero::PetersFiniteState{N,TF}, u, p, t) where {N,TF}
    # indices for extracting state variables
    iq = SVector{4}(1:4)
    iλ = SVector{N}(5:4+N)
    # separate aerodynamic and structural states
    q = u[iq]
    λ = u[iλ]
    # extract structural state variables
    h, θ, hdot, θdot = q
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # extract model constants
    bbar = aero.b
    # calculate induced flow velocity
    λ0 = 1/2 * bbar'*λ
    # calculate lift (excluding state rate terms)
    L = 2*pi*ρ*U*b*(hdot + U*θ + (b/2-a*b)*θdot) + pi*ρ*b^2*(U*θdot - λ0)
    # calculate moment (excluding state rate terms)
    M = -pi*ρ*b^3*U*θdot
    # return portion of inputs that is not dependent on the state rates
    return SVector(L, M, θdot, 0, 0)
end

nothing #hide
```

### Input Jacobian

The jacobian of ``f_y`` with respect to the state variables is calculated using the [`get_input_state_jacobian`](@ref) function for out-of-place combined input functions and [`get_input_state_jacobian!`](@ref) function for in-place combined inputs functions.

```@example developer
function get_Jrs(stru::TypicalSection, aero::PetersFiniteState, u, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip, a, b, U, ρ = p
    # extract model constants
    bbar = aero.b
    # create zero row vector with length nλ
    zλ = zero(aero.bbar')
    # compute jacobian sub-matrices
    Jrs = 2*pi*ρ*b*U*(@SMatrix [0 U 1 b-a*b; 0 0 0 -b^2/2])
    Jra = -pi/2*ρ*b^2*vcat(bbar', zero(bbar)')
    Jds = @SMatrix [0 0 0 1; 0 0 0 0; 0 0 0 0]
    Jda = vcat(zλ, zλ, zλ)
    # return jacobian
    return vcat(hcat(Jrs, Jra), hcat(Jds, Jda))
end
nothing #hide
```

## Model Ordering

For this package, we assume the following ordering for models, state variables, and parameters.

1. Rigid Body (when present)
2. Structural
3. Aerodynamics

We chose this ordering in order to place the most intuitive state variables first in the combined state variable vector.

## Avoiding Mass Matrices

In order to take advantage of as many features of the DifferentialEquations package as possible (including local sensitivity analysis) we recommend that the governing differential equations for models be reformulated to avoid using mass matrices whenever possible.
