# Theory

## Aerodynamic Models

### Two-Dimensional Models

#### [`Steady`](@ref)

This model is a steady linear 2D aerodynamic model derived from thin airfoil theory.  The equations for the lift and quarter-chord moment per unit span are:
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} \\
{M'}_{\frac{1}{4}} = 0
```
where ``a_0`` is the lift curve slope, ``\rho`` is the air density, ``u`` is the local freestream velocity in the chordwise direction, ``b`` is the semichord, and ``\alpha_\text{eff}`` is the effective angle of attack.  For this model, the effective angle of attack is
```math
\alpha_\text{eff} = -\frac{v}{u}
```
where ``v`` is the local freestream velocity in the airfoil normal direction and ``\alpha_0`` is the zero lift angle of attack.

For coupling with structural models, it is convenient to be able to define the lift and moment at an arbitrary reference location.  Defining the reference location to be ``a b`` aft of the semichord, the lift and moment may be expressed as
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} \\
M' = b \left(\frac{1}{2} + a \right) L'
```

#### [`QuasiSteady`](@ref)

This model is a quasi-steady linear 2D aerodynamic model derived from thin airfoil theory.  The equations for the lift and quarter-chord moment per unit span are:
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u \dot{\theta} - a b \ddot{\theta} \right) \\
{M'}_{\frac{1}{4}} = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right]
```
where ``a_0`` is the lift curve slope, ``\rho_\infty`` is the air density, ``u`` is the local freestream velocity in the chordwise direction, ``v`` is the local freestream velocity in the normal direction, ``\theta`` is the pitch angle, ``a`` defines the reference location, ``b`` is the semichord, and ``\alpha_\text{eff}`` is the effective angle of attack.  The reference location is located ``a b`` aft of the semichord.

The effective angle of attack ``\alpha`` for this model is given by
```math
\alpha_\text{eff} = -\frac{v}{u} + \frac{b}{u}\left( \frac{1}{2} - a \right) \dot{\theta} - \alpha_0
```
where ``\alpha_0`` is the zero lift angle of attack.

At the reference location, the lift and moment are
```math
L' = a_0 \rho_\infty u^2 b \alpha_\text{eff} + \pi \rho b^2 \left(-\dot{v} + u\dot{\theta} - a b \ddot{\theta} \right) \\
M' = -\pi \rho_\infty b^3 \left[ -\frac{1}{2}\dot{v} + u\dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right] + b \left(\frac{1}{2} + a \right) L'
```

#### [`Wagner`](@ref)

Wagner's function models the indicial response of aerodynamic loads under a sudden change in downwash ``w`` at the three-quarter's chord. The exact expression for Wagner's function is
```math
\phi(t) = \frac{2}{\pi} \int_0^\infty \frac{Re(C) \sin ( \omega (u/b) t  )}{\omega} d\omega
```
where ``C(\omega)`` is Theodorsen's function.  In many cases, approximate versions of Wagner's function are used rather than the exact expression,  of which one of the most common is the approximation of Wagner's function provided by R. T. Jones
```math
\phi(t) = 1 - C_1 e^{-\varepsilon_1 (u/b) t} - C_2 e^{-\varepsilon_2 (u/b) t}
```
where ``C_1 = 0.165``, ``C_2 = 0.335``, ``\varepsilon_1 = 0.455``, and ``\varepsilon_2 = 0.3``.

Wagner's function may be used to model arbitrary airfoil motion using Duhamel's integral.  We start by modeling the increment in circulatory lift ``d L_c(t)`` at time ``t`` due to an increment in downwash ``d w(t)`` at earlier time ``\tau`` as
```math
\frac{d L_c'(t)}{a_0 \rho_\infty u b} =  \phi(t - \tau) d w(\tau)
```
where ``\phi(t)`` is the impulse response function, which in this case is R. T. Jones' approximation of Wagner's function.  Superimposing all previous impulse responses using Duhamel's integral yields the following expression for the instantaneous circulatory lift.
```math
\frac{L_c'}{a_0 \rho_\infty u b} = \int_{-\infty}^t d w(\tau) \phi(t - \tau) d\tau = w(0) \phi(t) + \int_{0}^t  d w(\tau) \phi(t - \tau) d \tau
```
We can transform this equation using integration by parts, yielding
```math
\frac{L_c'}{a_0 \rho_\infty u b} = w(t) \phi(0) - \int_{0}^t w(\tau) d\phi(t - \tau) d\tau
```
The integral in this expression may be expressed as a function of the aerodynamic states ``\lambda_1`` and ``\lambda_2``.
```math
\lambda_1 = C_1 \varepsilon_1 \frac{u}{b} \int_0^t w(\tau) e^{-\varepsilon_1 (u/b) (t - \tau)} d \tau
```
```math
\lambda_2 = C_2 \varepsilon_2 \frac{u}{b} \int_0^t w(\tau) e^{-\varepsilon_2 (u/b) (t - \tau)} d \tau
```
```math
\lambda_1 + \lambda_2 = - \int_0^t w(\tau) d\phi(t-\tau) d\tau
```
The expression for circulatory lift then reduces to
```math
\frac{L_c'}{a_0 \rho_\infty u b} = w(t) \phi(0) + \lambda_1 + \lambda_2
```
where the downwash at the three quarter's chord is given by
```math
w(t) = -v + b \left( \frac{1}{2} - a \right) \dot{\theta} - u\alpha_0
```
and the aerodynamic states variables ``\lambda_1`` and ``\lambda_2`` are described by the ordinary differential equations
```math
\dot{\lambda_1} = -\varepsilon_1 \frac{u}{b} \lambda_1 + C_1 \varepsilon_1 \frac{u}{b} w(t) \\
\dot{\lambda_2} = -\varepsilon_2 \frac{u}{b} \lambda_2 + C_2 \varepsilon_2 \frac{u}{b} w(t)
```

The same lift and moment expressions are used as in the quasisteady model, but with the new effective angle of attack
```math
\alpha_\text{eff} = \left(-\frac{v}{u} + \frac{b}{u} \left( \frac{1}{2} - a \right) \dot{\theta} - \alpha_0 \right) \phi(0) + \frac{\lambda_1}{u} + \frac{\lambda_2}{u}
```

#### [`Peters`](@ref)

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

### Three-Dimensional Models

#### [`LiftingLine`](@ref)

Two-dimensional aerodynamic models may be applied in the context of a three-dimensional analysis by applying these models at multiple chordwise sections along the span of one or more lifting surfaces.  This type of model is applicable when spanwise flow effects are negligible, which is often the case for high aspect ratio wings.

The lifting line model implemented in this package assumes that the aerodynamics of each section is independent of the aerodynamics of the other sections, except as coupled through other models.  The state variables and inputs for this model correspond to the state variables, inputs, and parameters of each of the two-dimensional aerodynamic models, concatenated.  Rate equations are also concatenated.  

When coupled with a structural model, aircraft linear and angular accelerations are obtained from the structural model and transformed into the (deformed) local beam frame using an appropriate transformation matrix.  The local freestream velocities/accelerations and pitch rates/accelerations are then defined by a subset of the transformed linear and angular accelerations and cross-flow effects are neglected.   An inverse transformation may then be performed to transform the local aerodynamic forces/moments into the reference frame used by the structural model.

#### Vortex Lattice Method [`VLM`](@ref)

## Structural Models

### Two-Dimensional Models

#### [`TypicalSection`](@ref)

![](typical-section.svg)

The equations of motion for this model are
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M
```
where ``a`` is the normalized distance from the semichord to the reference point, ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_P`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M_\frac{1}{4}`` is the quarter-chord moment per unit span.

### Three-Dimensional Models

#### [`RigidBody`](@ref)

#### [`GEBT`](@ref)
