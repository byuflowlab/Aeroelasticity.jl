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
    MyTypicalSection

Typical section structural model with state variables ``q = \\begin{bmatrix} h &
θ & \\dot{h} & \\dot{\\theta} \\end{bmatrix}^T``, structural parameters ``p_s =
\\begin{bmatrix} a & b & k_h & k_\\theta & m & x_\\theta & I_P \\end{bmatrix}^T``,
and aerodynamic loads ``r = \\begin{bmatrix} L & M_\\frac{1}{4} \\end{bmatrix}^T``
"""
struct MyTypicalSection end

nothing #hide
```

Model parameters may be stored either as fields of the newly defined type or passed directly to the solver.  In this case, we choose to define all of the structural parameters as elements of the parameter vector.  In general, model constants should be stored as fields of the struct, whereas parameters that may change should be passed to the solver through the parameter vector.

Note that the state, parameter, and aerodynamic load variable identities should be documented in the docstring of the type definition since the type definition provides the primary (and sometimes only) source of documentation for a given model.

### Defining Model Properties

We now need to define a few model properties.

To indicate the number of state variables in our model, we define a new method for the [`number_of_states`](@ref) function.  In our case, we have four structural state variables ``q = \begin{bmatrix} h & \theta & \dot{h} & \dot{\theta} \end{bmatrix}^T ``.  This method must be defined for all models.

```@example developer
number_of_states(::Type{MyTypicalSection}) = 4
nothing #hide
```

To indicate the number of aerodynamic load inputs, we define a new method for the [`number_of_inputs`](@ref) function.  In our case, we have two aerodynamic loads ``r = \begin{bmatrix} L & M_\frac{1}{4} \end{bmatrix}^T``.  This method must be defined for all models.

```@example developer
number_of_inputs(::Type{MyTypicalSection}) = 2
nothing #hide
```

To indicate the number of parameters, we define a new method for the [`number_of_parameters`](@ref) function.  In our case, we have seven parameters ``p_s = \begin{bmatrix} a & b & k_h & k_\theta & m & x_\theta & I_P \end{bmatrix}^T``.  This method must be defined for all models.

```@example developer
number_of_parameters(::Type{MyTypicalSection}) = 7
nothing #hide
```

To indicate whether we plan to use in-place or out-of-place functions, we define a new method for the [`inplaceness`](@ref) function.  In general, for performance reasons, in-place functions are preferred.  The one exception is for models with small amounts of state variables, in which case the preferred approach is to use static arrays with out-of-place functions.  For this model, we use the latter approach.  This method must be defined for all models.

```@example developer
inplaceness(::Type{MyTypicalSection}) = OutOfPlace()
nothing #hide
```

To indicate mass matrix properties, we define a new method for the [`mass_matrix_type`](@ref) function.  This method must be defined for all models.

```@example developer
mass_matrix_type(::Type{MyTypicalSection}) = Varying()
nothing #hide
```

To indicate the properties of the jacobian of the state rates with respect to the state variables, we define a new method for the [`state_jacobian_type`](@ref) function.

```@example developer
state_jacobian_type(::Type{MyTypicalSection}) = Varying()
nothing #hide
```

To indicate the properties of the jacobian of the state rates with respect to the inputs, we define a new method for the [`input_jacobian_type`](@ref) function.

```@example developer
input_jacobian_type(::Type{MyTypicalSection}) = Varying()
nothing #hide
```

To indicate whether the state rates are linearly dependent on the aerodynamic loads, or in other words whether the governing structural equations can be expressed as  
```math
M_s(q,p_s,t)\dot{q} = f_s(q,p_s,t) + D_s(q,p_s,t)r
```
we define a new method for the [`input_dependence_type`](@ref) function.  Having a linear load dependence allows for greater flexibility when coupling the structural model with an aerodynamic model.

```@example developer
input_dependence_type(::Type{MyTypicalSection}) = Linear()
nothing #hide
```

### Mass Matrix Equation

For out-of-place models, the mass matrix is calculated using the [`get_mass_matrix`](@ref) function.  For in-place models the mass matrix is calculated using the [`get_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `q`, `r`, `p`, and `t` arguments.  

```@example developer
function get_mass_matrix(::MyTypicalSection, q, r, p, t)
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # calculate mass matrix
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
end
nothing #hide
```

### State Rate Equation

The right hand side of the governing structural differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models and [`get_rates!`](@ref) function for in-place models.  

```@example developer
function get_rates(::MyTypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # calculate state rates
    return SVector(hdot, θdot, -kh*h - L, -kθ*θ + M + (b/2+a*b)*L)
end

nothing #hide
```

### State Rate Jacobian

The jacobian of the right hand side of the governing structural equations with respect to the state variables is calculated using the [`get_state_jacobian`](@ref) function for out-of-place models and [`get_state_jacobian!`](@ref) function for in-place models.

```@example developer
function get_state_jacobian(::MyTypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

nothing #hide
```

### Input Jacobian

The jacobian of the right hand side of the governing structural equations with respect to the inputs is defined using the [`get_input_jacobian`](@ref) function.  There is no out-of-place form for this function, however, it may be constructed as a either a linear map (if large) or static array (if small) in order to avoid allocations.  

```@example developer
function get_input_jacobian(::MyTypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0; 0 0; -1 0; b/2+a*b 1]
end

nothing #hide
```

### Structural Model Code

Putting it all together, a complete representation of our typical section model for use with this package may be defined using the following block of code.

```@example developer-combined
"""
    MyTypicalSection <: AbstractModel

Typical section structural model with state variables ``q = \\begin{bmatrix} h &
θ & \\dot{h} & \\dot{\\theta} \\end{bmatrix}^T``, structural parameters ``p_s =
\\begin{bmatrix} a & b & k_h & k_\\theta & m & x_\\theta & I_P \\end{bmatrix}^T``,
and aerodynamic loads ``r = \\begin{bmatrix} L & M_\\frac{1}{4} \\end{bmatrix}^T``
"""
struct MyTypicalSection <: AbstractModel end

# --- Traits --- #

number_of_states(::Type{MyTypicalSection}) = 4
number_of_inputs(::Type{MyTypicalSection}) = 2
number_of_parameters(::Type{MyTypicalSection}) = 7
inplaceness(::Type{MyTypicalSection}) = OutOfPlace()
mass_matrix_type(::Type{MyTypicalSection}) = Varying()
state_jacobian_type(::Type{MyTypicalSection}) = Varying()
input_jacobian_type(::Type{MyTypicalSection}) = Varying()
input_dependence_type(::Type{MyTypicalSection}) = Linear()

# --- Methods --- #

function get_mass_matrix(::MyTypicalSection, q, r, p, t)
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # calculate mass matrix
    return @SMatrix [1 0 0 0; 0 1 0 0; 0 0 m m*b*xθ; 0 0 m*b*xθ Ip]
end

function get_rates(::MyTypicalSection, q, r, p, t)
    # extract structural states
    h, θ, hdot, θdot = q
    # extract aerodynamic loads
    L, M = r
    # extract structural parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # calculate state rates
    return SVector(hdot, θdot, -kh*h - L, -kθ*θ + M + (b/2+a*b)*L)
end

function get_state_jacobian(::MyTypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0 1 0; 0 0 0 1; -kh 0 0 0; 0 -kθ 0 0]
end

function get_input_jacobian(::MyTypicalSection, q, r, p, t)
    # extract parameters
    a, b, kh, kθ, m, xθ, Ip = p
    # return jacobian
    return @SMatrix [0 0; 0 0; -1 0; b/2+a*b 1]
end

nothing #hide
```

## Defining an Aerodynamic Model

For demonstrating how to create a new aerodynamic model, we use Peters' finite state aerodynamic model.

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

### Theory

Aerodynamic state variables in this package satisfy the ordinary differential equation (or differential algebraic equation in mass matrix form)
```math
M_{a}(λ,d,p_a,t)\dot{\lambda} = f_a(λ,d,p_a,t)
```
where ``M_a`` is a function which defines the aerodynamic mass matrix, ``f_a`` is a function which defines the mass matrix multiplied aerodynamic state rates, ``\lambda`` is a vector of structural states, ``d`` is a vector of structural deflections, ``p_a`` is a vector of aerodynamic parameters, and ``t`` is the current time.

The equations of motion for Peters' finite state model, when expressed in the form expected by this package are:
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
    MyPeters{N,TF,SV,SA} <: AbstractModel

Peters' finite state model with `N` state variables, inputs ``d = \\begin{bmatrix}
\\dot{\\theta} & \\ddot{h} & \\ddot{\\theta}\\end{bmatrix}^T`` and parameters
``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``
"""
struct MyPeters{N,TF,TV<:SVector{N,TF},TA<:SMatrix{N,N,TF}} <: AbstractModel
    A::TA
    b::TV
    c::TV
end

nothing #hide
```

Here `N` is the number of aerodynamic state variables and `TF` is the floating point type used to represent the constant matrices/vectors ``\bar{A}``, ``\bar{b}``, and ``\bar{c}``.

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

To indicate the number of inputs, we define a new method for the [`number_of_inputs`](@ref) function.  In our case, we have three inputs ``d = \begin{bmatrix} \dot{\theta} & \ddot{h} & \ddot{\theta} \end{bmatrix}^T``.  This method must be defined for all models.

```@example developer
number_of_inputs(::Type{<:MyPeters}) = 3
nothing #hide
```

To indicate the number of parameters, we define a new method for the [`number_of_parameters`](@ref) function.  In our case, we have four parameters ``p_a = \begin{bmatrix} a & b & U & \rho \end{bmatrix}^T``.  This method must be defined for all models.

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
state_jacobian_type(::Type{<:MyPeters}) = Varying()
nothing #hide
```

To indicate the properties of the jacobian of the state rates with respect to the inputs, we define a new method for the [`input_jacobian_type`](@ref) function.

```@example developer
input_jacobian_type(::Type{<:MyPeters}) = Varying()
nothing #hide
```

To indicate whether the state rates are linearly dependent on the inputs, or in other words whether the governing aerodynamic equations can be expressed as  
```math
M_a(\lambda,p_a,t)\dot{\lambda} = f_a(\lambda,p_s,t) + D_a(\lambda,p_s,t)d
```
we define a new method for the [`input_dependence_type`](@ref) function.

```@example developer
input_dependence_type(::Type{<:MyPeters}) = Linear()
nothing #hide
```

### Mass Matrix Equation

For out-of-place models, the mass matrix is calculated using the [`get_mass_matrix`](@ref) function.  For in-place models the mass matrix is calculated using the [`get_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `λ`, `d`, `p`, and `t` arguments.  

```@example developer
get_mass_matrix(model::MyPeters) = model.A
nothing #hide
```

### State Rate Equation

The right hand side of the governing aerodynamic differential equations is calculated using the [`get_rates`](@ref) function for out-of-place models and [`get_rates!`](@ref) function for in-place models.  

```@example developer
function get_rates(model::MyPeters{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract structural deflections
    θdot, hddot, θddot = d
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return cbar*(hddot + U*θdot + (b/2-a*b)*θddot) - U/b*λ
end

nothing #hide
```

### State Rate Jacobian

The jacobian of the right hand side of the governing aerodynamic equations with respect to the state variables is calculated using the [`get_state_jacobian`](@ref) function for out-of-place models and [`get_state_jacobian!`](@ref) function for in-place models.

```@example developer
function get_state_jacobian(model::MyPeters, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return -U/b*Diagonal(one.(cbar))
end

nothing #hide
```

### Structural Deflection Jacobian

The jacobian of the right hand side of the governing aerodynamic equations with respect to the structural deflections is defined using the [`get_input_jacobian`](@ref) function.

```@example developer
function get_input_jacobian(model::MyPeters, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return hcat(U*cbar, cbar, (b/2-a*b)*cbar)
end

nothing #hide
```

### Aerodynamic Model Code

Putting it all together, a complete representation Peters' finite state aerodynamic model may be defined using the following block of code.

```@example developer-combined
"""
    MyPeters{N,TF,SV,SA} <: AbstractModel

Peter's finite state model with `N` state variables, inputs ``d = \\begin{bmatrix}
\\dot{\\theta} & \\ddot{h} & \\ddot{\\theta}\\end{bmatrix}^T`` and parameters
``p_a = \\begin{bmatrix} a & b & U & \\rho \\end{bmatrix}^T``
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
number_of_inputs(::Type{<:MyPeters}) = 3
number_of_parameters(::Type{<:MyPeters}) = 4
inplaceness(::Type{<:MyPeters}) = OutOfPlace()
mass_matrix_type(::Type{<:MyPeters}) = Constant()
state_jacobian_type(::Type{<:MyPeters}) = Varying()
input_jacobian_type(::Type{<:MyPeters}) = Varying()
input_dependence_type(::Type{<:MyPeters}) = Linear()

# --- Methods --- #

get_mass_matrix(model::MyPeters) = model.A

function get_rates(model::MyPeters{N,TF,SV,SA}, λ, d, p, t) where {N,TF,SV,SA}
    # extract aerodynamic states as statically sized vector
    λ = SVector{N}(λ)
    # extract structural deflections
    θdot, hddot, θddot = d
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # calculate rates
    return cbar*(hddot + U*θdot + (b/2-a*b)*θddot) - U/b*λ
end

function get_state_jacobian(model::MyPeters, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # jacobian with respect to aerodynamic states
    return -U/b*Diagonal(one.(cbar))
end

function get_input_jacobian(model::MyPeters, λ, d, p, t)
    # extract parameters
    a, b, U, ρ = p
    # extract model constants
    cbar = model.c
    # return jacobian
    return hcat(U*cbar, cbar, (b/2-a*b)*cbar)
end

nothing #hide
```

## Coupling Aerodynamic and Structural Models

For demonstrating how to couple aerodynamic and structural models together, we use the typical section structural model and Peters' finite state aerodynamic model.  To couple these two models together, we need to define the model inputs (aerodynamic loads and structural deflections) as functions of the states, parameters, and time.

### Theory

We assume the aerodynamic loads and structural deflections may be expressed as a function of the aerodynamic and structural state variables and parameters, as well as the current time.

```math
d = f_d(u,p,t) \quad r = f_r(u,p,t) \\
```
where
```math
u = \begin{bmatrix} \lambda & q \end{bmatrix}^T \quad
p = \begin{bmatrix} p_a & p_s \end{bmatrix}^T \quad
```

If the state rates of a given structural model are linearly dependent on the aerodynamic loads, we can expand the expression which defines the aerodynamic loads to
```math
r = f_r(u,p,t) - M_{r s}(u,p,t) \dot{q} -
M_{r a}(u,p,t) \dot{\lambda}
```
where ``M_{r s}`` is a function which defines the (negative) jacobian of the aerodynamic loads with respect to the structural state rates and ``M_{r a}`` is a function which defines the (negative) jacobian of the aerodynamic loads with respect to the aerodynamic state rates.  ``f_r`` is a function which defines the portion of the aerodynamic loads which is independent of the structural and aerodynamic state rates.

If the state rates of a given aerodynamic model are linearly dependent on the structural deflections, we can expand the expression which defines the structural deflections to
```math
d = f_d(u,p,t) - M_{d s}(u,p,t) \dot{q} -
M_{d a}(u,p,t) \dot{\lambda}
```
where ``M_{d s}`` is a function which defines the (negative) jacobian of the structural deflections with respect to the structural state rates and ``M_{d a}`` is a function which defines the (negative) jacobian of the structural deflections with respect to the aerodynamic state rates.  ``f_d`` is a function which defines the portion of the structural deflections which is independent of the structural and aerodynamic state rates.

In the most general case, the coupled system of equations may be defined as
```math
M_u(u,p,t) \dot{u} = f_u(u,p,t)
```
where
```math
f_u = \begin{bmatrix} f_a \\f_s \end{bmatrix} \quad
M_{u} = \begin{bmatrix} M_a & 0 \\ 0 & M_s \end{bmatrix} + \begin{bmatrix} D_a & 0 \\ 0 & D_s \end{bmatrix} \begin{bmatrix} M_{da} & M_{ds} \\ M_{ra} & M_{rs} \end{bmatrix} \quad
```
The associated jacobian is
```math
\frac{\partial f_u}{\partial u} = \begin{bmatrix} \frac{\partial f_a}{\partial \lambda} \end{bmatrix} & 0 \\ 0 & \frac{\partial f_s}{\partial q} + \begin{bmatrix} D_a & 0 \\ 0 & D_s \end{bmatrix} \begin{bmatrix} \frac{\partial f_d}{\partial \lambda} & \frac{\partial f_d}{\partial q} \\ \frac{\partial f_r}{\partial \lambda} & \frac{\partial f_r}{\partial q} \end{bmatrix}
```

If we introduce the combined input function
```math
y = f_y(u, p, t) - M_y(u, p, t) \dot{u}
```
where
```math
y = \begin{bmatrix} d \\ r \end{bmatrix} \quad f_y = \begin{bmatrix} f_d \\ f_r \end{bmatrix} \quad M_y = \begin{bmatrix} M_{da} & M_{ds} \\ M_{ra} & M_{rs} \end{bmatrix}
```
the mass matrix and jacobian expressions may be shortened to
```math
M_{u} = \begin{bmatrix} M_a & 0 \\ 0 & M_s \end{bmatrix} + \begin{bmatrix} D_a & 0 \\ 0 & D_s \end{bmatrix} M_y
```
```math
\frac{\partial f_u}{\partial u} =
\begin{bmatrix}
\frac{\partial f_a}{\partial \lambda} & 0 \\
0 & \frac{\partial f_s}{\partial q}
\end{bmatrix} +
\begin{bmatrix}
D_a & 0 \\
0 & D_s
\end{bmatrix} \frac{\partial f_y}{\partial u}
```

We have already implemented all of the functions in these expressions with the exception of those associated with the combined input function, so all that is required to couple the aerodynamic and structural models together is to define functions associated with the combined input function.

#### Aerodynamic Model Inputs

The inputs expected by Peters' finite state model are ``d = \begin{bmatrix} \dot{\theta} & \ddot{h} & \ddot{\theta} \end{bmatrix}^T``.  These structural deflections correspond to a subset of the structural states and corresponding rates.

The structural deflections may be expressed in the form expected by this package as
```math
d = J_{ds} q - M_{ds} \dot{q}
```
where
```math
J_{ds} = \begin{bmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \quad
M_{ds} = \begin{bmatrix} 0 & 0 & 0 & 0 \\ 0 & 0 & -1 & 0 \\ 0 & 0 & 0 & -1 \end{bmatrix}
```

#### Structural Model Inputs

The aerodynamic loads expected by the typical section model are the lift and quarter-chord moment.  The lift and quarter-chord moment, as calculated using Peters' finite state model are
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

For out-of-place input functions, the mass matrix is calculated using the [`get_input_mass_matrix`](@ref) function.  For in-place combined input functions, the mass matrix is calculated using the [`get_input_mass_matrix!`](@ref) function.  For constant mass matrices, these functions are called without the `u`, `p`, and `t` arguments.

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

The jacobian of ``f_y`` with respect to the state variables is calculated using the [`get_input_state_jacobian`](@ref) function for out-of-place combined input functions and [`get_input_state_jacobian!`](@ref) function for in-place combined inputs functions.

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
