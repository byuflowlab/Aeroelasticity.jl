# Models

## Aerodynamic Models

### Two-Dimensional Models

#### Steady [`Steady`](@ref)

This model is a steady linear 2D aerodynamic model derived from thin airfoil theory.  The equations for the lift and quarter-chord moment per unit span are:
```math
L' = a_0*\rho_\infty*U^2*b*(\alpha-\alpha_0) \\
{M'}_{\frac{1}{4}} = 0
```
where ``a_0`` is the lift curve slope, ``\rho`` is the air density, ``U`` is the freestream velocity, ``b`` is the semichord, ``\alpha`` is the effective angle of attack, and ``\alpha_0`` is the zero lift angle of attack.  For this model, the effective angle of attack ``\alpha`` is assumed to be equal to the pitch angle ``\theta``.

For coupling with structural models, it is convenient to be able to define the lift and moment at an arbitrary reference location.  Defining the reference location to be ``a b`` aft of the semichord, the lift and moment may be expressed as
```math
L' = a_0*\rho_\infty*U^2*b*(\alpha-\alpha_0) \\
M' = b \left(\frac{1}{2} + a \right) L'
```

#### Quasi-Steady [`QuasiSteady`](@ref)

This model is a quasi-steady linear 2D aerodynamic model derived from thin airfoil theory.  The equations for the lift and quarter-chord moment per unit span are:
```math
L' = \a_0*\rho_\infty*U^2*b*(\alpha - \alpha_0) + \pi \rho b^2 \left(\ddot{h} + U\dot{\theta} - a b \ddot{\theta}) \\
{M'}_{\frac{1}{4}} = -\pi \rho_\infty b^3 \left[ \frac{1}{2}\ddot{h} + U\dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right]
```
where ``a_0`` is the lift curve slope, ``\rho`` is the air density, ``U`` is the freestream velocity, ``a`` defines the reference location, ``b`` is the semichord, ``h`` is the plunge, ``\alpha`` is the effective angle of attack, ``\alpha_0`` is the zero lift angle of attack, and ``\theta`` is the pitch angle.  The reference location is located ``a b`` aft of the semichord.  The lift contains both circulatory and noncirculatory terms, whereas the pitching moment about the quarter chord is entirely noncirculatory.

The effective angle of attack for this model is
```math
\alpha = \theta + \frac{\dot{h}}{U} + \frac{b}{U}\left( \frac{1}{2} - a \right) \dot{\theta}
```

At the reference location, the lift and moment are
```math
L' = \a_0*\rho_\infty*U^2*b*(\alpha - \alpha_0) + \pi \rho b^2 \left(\ddot{h} + U\dot{\theta} - a b \ddot{\theta}) \\
{M'}_{\frac{1}{4}} = -\pi \rho_\infty b^3 \left[ \frac{1}{2}\ddot{h} + U\dot{\theta} + b \left( \frac{1}{8} - \frac{a}{2} \right) \ddot{\theta} \right] + b \left(\frac{1}{2} + a \right) L'
```

#### Wagner's Function [`Wagner`](@ref)

We can model the increment in circulatory lift ``\delta L_c(t)`` at time ``t`` due to an increment in effective angle of attack ``\delta \alpha`` at time ``\tau`` as
```math
\frac{\delta L_c'}{a_0 \rho_\infty U^2 b} =  \phi(t - \tau) \delta \alpha(\tau)
```
where ``\phi(t)`` is the impulse response function. Superimposing all previous impulse responses using Duhamel's integral yields the following expression for the instantaneous circulatory lift.
```
\frac{L_c'}{a_0 \rho_\infty U^2 b} = \int_{-\infty}^t \phi(t - \tau) \delta \alpha(\tau)
```

Wagner's function is given by
```math
\phi(t) = 1/2 - \int_0^\infty \frac{Im\left[ e^{iwt} C(\omega) \right]}{\omega} d\omega
```
where ``C(\omega)`` is Theodorsen's function. This function models the indicial response of aerodynamic loads under a sudden change in downwash ``w = U \sin(\alpha) \approx U \alpha`` and is often approximated as:
```math
ϕ(t) = 1 - 0.165 e^-0.455t - 0.335 e^-0.3t
```



Wagner's

In order to model the delay between changes in the downwash of an airfoil ``\delta w(t)`` and corresponding changes in circulatory lift ``\delta L_c(t)'``, the  may be introduced.




To model the time history of the circulatory lift for arbitrary time-varying motion, we can approximate the time history of the downwash as a sum of several step signals and use Duhamel's integral.
```
L_c' = \int_{-\infty}^t \phi(t - \tau) \delta L_c'(0)
```

Duhamel's intergral may then be used to approximate the time history of the downwash as a sum of several step signals and use Duhamel's integral


In the quasi-steady aerodynamic model, we assume that instantaneous angle of attack perturbations correspond to instantenous changes in the airfoil lift.
 correspond to instantaneous changes in the
The quasi-steady aerodynamic model assumes that lift is solely a  



We start by introduc

Wagner's function models the indicial response of aerodynamic loads under a sudden change in downwash.  

To

```math
dw = dw_{qs}(t)\phi(t - \tau)
```

To model the time history of the downwash for arbitrary time-varying motion, we can approximate the time history of the downwash as a sum of several step signals and use Duhamel's integral with Wagner's function.  The result is the following expression for the downwash at the three-quarter chord point:

```math
w(t) = w(0)ϕ(t) + \int_0^\partial{w}\partial{\tau}
```

 Using  To use Wagner's function in the context of an arbitrary time-varying model, we model the time history of the downwash as

Circulatory lift calculated using Wagner's function is
```math
L = \a_0*\rho_\infty*U*b*(\phi(t)
```
where ``\phi(t)`` is Wagner's function.



Using Wagner's function, the circulat

To use Wagner's function, the circulatory lift in  is defined as

To use Wagner's function in the context of a state space model, the effective angle of attack is modified to account for the
```math
\alpha = \theta + \frac{\dot{h}}{U} + \frac{b}{U}\left( \frac{1}{2} - a \right) \dot{\theta}
```


Defining the lift and quarter chord moment as:


#### Peter's Finite State Model [`Peters`](@ref)

### Three-Dimensional Models

#### Lifting Line [`LiftingLine`](@ref)

#### Vortex Lattice Method [`VLM`](@ref)

## Structural Models

### Two-Dimensional Models

#### Typical Section [`TypicalSection`](@ref)

![](typical-section.svg)

The equations of motion for this model are
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M_{\frac{1}{4}} + b \left( \frac{1}{2} + a \right) L
```
where ``a`` is the normalized distance from the semichord to the reference point, ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_P`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M_\frac{1}{4}`` is the quarter-chord moment per unit span.

### Three-Dimensional Models

#### Rigid Body [`RigidBody`](@ref)

#### Geometrically Exact Beam Theory [`GEBT`](@ref)
