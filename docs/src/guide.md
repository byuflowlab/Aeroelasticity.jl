# [Getting Started](@id guide)

The governing equations for many unsteady systems may be described by the first-order
implicit differential equation
```math
   0 = \bm{f}(\bm{\dot{x}}, \bm{x}, \bm{y}, \bm{p}, t)
```
where ``\bm{f}`` is a vector-valued residual function, ``\dot{()}`` denotes the time derivative,
``\bm{x}`` is a vector of state variables, ``\bm{y}`` is a vector of time-varying inputs,
``\bm{p}`` is a vector of time-invariant parameters, and ``t`` is the current time. If we
concatenate the governing equations associated with any number of unsteady models described
by this equation together, we obtain the concatenated model
```math
\bm{0} = \bm{\mathcal{F}}(\bm{\dot{X}}, \bm{X}, \bm{Y}, \bm{P}, t)
```
where
```math
\bm{\mathcal{F}} = \begin{bmatrix}
\bm{f_1} \\
\bm{f_2} \\
\vdots \\
\bm{f_n}
\end{bmatrix} \quad
\bm{X} = \begin{bmatrix}
\bm{x_1} \\
\bm{x_2} \\
\vdots \\
\bm{x_n}
\end{bmatrix} \quad
\bm{Y} = \begin{bmatrix}
\bm{y_1} \\
\bm{y_2} \\
\vdots \\
\bm{y_n}
\end{bmatrix} \quad
\bm{P} = \begin{bmatrix}
\bm{p_1} \\
\bm{p_2} \\
\vdots \\
\bm{p_n}
\end{bmatrix}
```
The concatenated model is a decoupled model, since the states, inputs, and parameters of
each submodel are not influenced in any way by the states, inputs, and parameters of the
other submodels. To couple these models, we introduce the following coupling function,
which defines the inputs of the concatenated model as functions of its rates, states, and
parameters, as well as the current time.
```math
\bm{y} = \bm{\mathcal{G}}(\bm{\dot{X}}, \bm{X}, \bm{P}, t)
```
Defining this coupling function allows us to define the governing equations for a general
monolithic coupled model.
```math
\bm{\tilde{\mathcal{F}}}(\bm{\dot{X}}, \bm{X}, \bm{P}, t) = \bm{\mathcal{F}}(\bm{\dot{X}}, \bm{X}, \bm{\mathcal{G}}(\bm{\dot{X}}, \bm{X}, \bm{P}, t), \bm{P}, t) = \bm{0}
```

## Implementation

Our coupling methodology is implemented by the following `CoupledModel` constructor.

```@docs
CoupledModel(submodels, parameters)
```

In addition to constructing a coupled model, this constructor automatically detects the
sparsity of the resulting rate and state jacobian matrices.  This information is then used
when evaluating the following jacobian evaluation functions.

```@docs
rate_jacobian!
state_jacobian!
```

The residual may be evaluated using the following function

```@docs
residual!
```

Alternatively, the residual may be evaluated by calling the coupled model as a function.

```@docx
(coupled_model::CoupledModel)(resid, dx, x, p, t)
```

## Analyses

Once a coupled model has been constructed, the system may be linearized using the
[`linearize`](@ref) function.

```@docs
linearize
```

Eigenvalues and eigenvectors may then be computed using the `dense_eigen` or
`sparse_eigen` functions.

```@docs
dense_eigen
sparse_eigen
```

Correlating eigenmodes with modes from previous iterations may be done using the
`correlate_eigenmodes` function

```@docs
correlate_eigenmodes
```

[`DifferentialEquations`](https://github.com/SciML/DifferentialEquations.jl) may be used to 
find a steady state or time marching solution.  To facilitate these analyses, this package 
provides the following specialized constructors for `ODEFunction` and `DAEFunction`.

```@docs
Aeroelasticity.SciMLBase.ODEFunction(model::CoupledModel, p0)
Aeroelasticity.SciMLBase.DAEFunction(model::CoupledModel)
```

## Aeroelasticity

The primary intended purpose of the coupling methodology presented in this package 
is to faciltate constructing aeroelastic systems.  To this end, a number of aerodynamic
and structural models have been implemented in this package which may be used to construct
various aeroelastic models.  Each model is defined as a callable struct, with a calling
signature that matches the format expected by the [`CoupledModel`](@ref) constructor. For
more details about the theoretical background behind each model, we refer readers to the 
relevant pages of the documentation.  For more details about how to use these models to 
construct aeroelastic systems, see the examples.
