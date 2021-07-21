# Peters' Finite State Model

## Theory

For Peter's finite state model, an additional term is added to the expression for the effective angle of attack from the quasi-steady model to account for induced velocity.
```math
\alpha = - \frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \dot{\theta} + \frac{\lambda_0}{u}
```

The induced velocity ``\lambda_0`` is approximated from a set of N induced-flow states ``\lambda_1, \lambda_2, \dots, \lambda_N`` as
```math
\lambda \approx \frac{1}{2} \sum_{n=1}^N b_n \lambda_n
```
The set of N first-order ordinary differential equations which govern the N finite aerodynamic states are derived by Peters as
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

The same lift and moment expressions are used as in the quasisteady model, but with the new effective angle of attack.

## Documentation

```@docs
Peters()
```
