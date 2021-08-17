# [Peters' Finite State Model](@id peters-finite-state-model)

## Theory

For Peter's finite state model, an additional term is added to the expression for the effective angle of attack from the [quasi-steady model](@ref quasi-steady-thin-airfoil-theory) to account for induced velocity.
```math
\alpha = - \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \omega + \frac{\lambda_0}{u} - \alpha_0
```

The induced velocity ``\lambda_0`` is approximated from a set of N induced-flow states ``\lambda_1, \lambda_2, \dots, \lambda_N`` as
```math
\lambda_0 \approx \frac{1}{2} \sum_{n=1}^N b_n \lambda_n
```
The set of N first-order ordinary differential equations which govern the N finite aerodynamic states are derived by Peters et al. as
```math
\bar{A} \lambda + \frac{u}{b} \lambda = \bar{c} \left[ -\dot{v} + u\omega + b \left(\frac{1}{2} - a \right) \dot{\omega} \right]
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

The same force and moment expressions are used as in the [quasi-steady model](@ref quasi-steady-thin-airfoil-theory), but with the new effective angle of attack provided above.

Note that while Peter's finite state model can accommodate any number of aerodynamic state variables, generally 3-10 state variables are used.

## Type Definition

```@docs
Peters
```

## Constructors

```@docs
Peters()
```

## Example Usage
 - [Aeroelastic Analysis of a Typical Section](@ref)
