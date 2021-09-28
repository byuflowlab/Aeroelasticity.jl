# Examples

## Aeroelastic Analysis of a Typical Section

```@setup typical-section
using Plots
pyplot()
nothing #hide
```

In this example, we demonstrate how to perform a two-dimensional aeroelastic analysis using a typical section model with two degrees of freedom.

![](typical-section.svg)

The equations of motion for this model are
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M
```
where ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_θ`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M`` is the moment per unit span about the reference point.

This package has a number of different two-dimensional aerodynamic models which we may use to model the aerodynamics of the typical section model.  These models include a steady-state model, a quasi-steady model, an unsteady aerodynamic model based on Wagner's function, and Peters' finite state aerodynamic model.  We will perform an aeroelastic analysis using each of these models and compare the results.

The non-dimensional parameters we will use are
```math
a = -1/5 \quad e = -1/10  \\
r^2 = \frac{I_P}{m b^2} \quad \sigma = \frac{\omega_h}{\omega_\theta} \\
\mu = \frac{m}{\rho_\infty \pi b^2} \quad V = \frac{U}{b \omega_\theta}
```
where ``a`` is the normalized distance from the semichord to the reference point, ``e`` is the normalized distance from the semichord to the center of mass, and ``\omega_h`` and ``\omega_\theta`` are the uncoupled natural frequencies, defined as follows.
```math
\omega_h = \sqrt{\frac{k_h}{m}} \quad \omega_\theta = \sqrt{\frac{k_\theta}{I_P}}
```

```@example typical-section-stability
using AerostructuralDynamics, DifferentialEquations, LinearAlgebra

# reduced velocity range
V = range(1e-6, 3.1, length=5000) # = U/(b*ωθ) (reduced velocity)

# non-dimensional parameters
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
cd0 = 0
cm0 = 0

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
aerodynamic_models = (Steady(), QuasiSteady(), Wagner(), Peters{6}())

# structural model
structural_model = TypicalSection()

# eigenvalue/eigenvector storage
λ = Vector{Matrix{ComplexF64}}(undef, length(aerodynamic_models))
Uλ = Vector{Array{ComplexF64,3}}(undef, length(aerodynamic_models))
Vλ = Vector{Array{ComplexF64,3}}(undef, length(aerodynamic_models))

# loop through each aerodynamic model
for (ia, aerodynamic_model) in enumerate(aerodynamic_models)

    # coupled model
    model = couple_models(aerodynamic_model, structural_model)

    # ode representation of the coupled model
    f = get_ode(model)

    # equilibrium rates are zero
    dx = zeros(number_of_states(model))

    # initial guess for equilibrium states
    x0 = zeros(number_of_states(model))

    # current time
    t = 0.0

    # eigenvalue/eigenvector storage
    nλ = number_of_states(model)
    λ[ia] = zeros(ComplexF64, nλ, length(V))
    Uλ[ia] = zeros(ComplexF64, nλ, nλ, length(V))
    Vλ[ia] = zeros(ComplexF64, nλ, nλ, length(V))

    # loop through each reduced frequency
    for i = 1:length(V)

        # set parameters
        p_aero = [a, b, a0, α0, cd0, cm0]
        p_stru = [kh, kθ, m, Sθ, Iθ]
        p_input = [U[i], ρ]
        p = vcat(p_aero, p_stru, p_input)

        # find equilibrium point
        x = solve(SteadyStateProblem(f, x0, p))

        # find corresponding inputs
        y = get_coupling_inputs(model, dx, x, p, t)

        # perform linear stability analysis
        λi, Uλi, Vλi = get_eigen(model, dx, x, y, p, t)

        # correlate eigenvalues
        if i > 1
            # previous left eigenvector matrix
            Uλpi = Uλ[ia][:,:,i-1]

            # current mass matrix
            Mi = get_rate_jacobian(model, dx, x, y, p, t)

            # use correlation matrix to correlate eigenmodes
            perm, corruption = correlate_eigenmodes(Uλpi, Mi, Vλi)

            # re-arrange eigenmodes
            λi = λi[perm]
            Uλi = Uλi[perm,:]
            Vλi = Vλi[:,perm]
        end

        # save eigenvalues/eigenvectors
        λ[ia][:,i] = λi
        Uλ[ia][:,:,i] = Uλi
        Vλ[ia][:,:,i] = Vλi
    end
end

nothing #hide
```

We now plot the results for each aerodynamic model.

```@example typical-section
using Plots
pyplot()

sp1 = plot(
    title = "Non-Dimensional Frequency",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    framestyle = :zerolines,
    titlefontsize = 14,
    guidefontsize = 14,
    legendfontsize = 11,
    tickfontsize = 11,
    legend = :topright,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true
    )

sp2 = plot(
    title = "Non-Dimensional Damping",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (-0.7, 0.605),
    ytick = -0.6:0.2:0.6,
    ylabel = "\$ \\frac{Γ}{\\omega_\\theta} \$",
    framestyle = :zerolines,
    titlefontsize = 14,
    guidefontsize = 14,
    legendfontsize = 11,
    tickfontsize = 11,
    legend = :topleft,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid = true
    )

labels = ["Steady", "Quasi-Steady", "Wagner", "Peters (N=6)"]

for ia = 1:length(aerodynamic_models)

    plot!(sp1, V, imag.(λ[ia][1,:])/ωθ,
        label = labels[ia],
        color = ia,
        markersize = 1,
        markerstrokewidth = 0,
        )

    for i = 2:size(λ[ia], 1)
        plot!(sp1, V, imag.(λ[ia][i,:])/ωθ,
            label = "",
            color = ia,
            markersize = 1,
            markerstrokewidth = 0,
            )
    end

    plot!(sp2, V, real.(λ[ia][1,:])/ωθ,
        label = labels[ia],
        color = ia,
        markersize = 1,
        markerstrokewidth = 0,
        )

    for i = 2:size(λ[ia], 1)
        plot!(sp2, V, real.(λ[ia][i,:])/ωθ,
            label = "",
            color = ia,
            markersize = 1,
            markerstrokewidth = 0,
            )
    end
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800))

savefig(p1, "typical-section-stability.svg") #hide

nothing #hide
```

![](typical-section-stability.svg)

Using the Wagner and/or Peters aerodynamic models yields a flutter reduced velocity around 2.2, while the steady and/or quasi-steady aerodynamic models predict significantly lower flutter velocities.  The aerodynamic state variables of the Wagner and Peters models allows these models to capture the impact of vortex shedding on the lift and drag of the profile, therefore we can expect these models to yield more accurate results than the steady-state and quasi-steady models.

We can visualize and/or create animations of the flutter mode with the help of custom plot recipes.

```@example typical-section

# flutter reduced velocity
Vf = 2.2

# aerodynamic model
ia = 4 # Peter's finite state model

# coupled model
model = couple_models(aerodynamic_models[ia], structural_model)

# flutter mode identity
it = argmin(abs.(V .- Vf))
iλ = argmin(abs.(real.(λ[ia][:,it])))

# flutter mode eigenvalue and eigenvector
λf = λ[ia][iλ,it]
vf = Vλ[ia][:,iλ,it]

# flutter mode state rates
dxf = zeros(number_of_states(model))

# flutter mode state variables
xf = zeros(number_of_states(model))

# flutter mode parameters
p_aero = [a, b, a0, α0, cd0, cm0]
p_stru = [kh, kθ, m, Sθ, Iθ]
p_input = [Vf*b*ωθ, ρ]
pf = vcat(p_aero, p_stru, p_input)

# animation time
t1 = -pi/abs(imag(λf))
t2 = pi/abs(imag(λf))

# eigenvector scaling
scaling = 0.5

# create animation
anim = @animate for t in range(t1, t2, length=50)

    xi = xf + scaling*real.(vf*exp(λf*t))

    plot(model, dxf, xi, pf, t)

end

# save animation
gif(anim, "typical-section-flutter-mode.gif")

nothing #hide
```

The non-dimensional parameters we use for this example match those used by Hodges and Pierce in "Introduction to Structural Dynamics and Aeroelasticity".  Hodges and Pierce performed the analysis using a steady state model and Peter's finite state model with six state variables.   The results presented here for the steady-state and Peters' finite state models match the results presented by Hodges and Pierce in "Introduction to Structural Dynamics and Aeroelasticity", which validates our implementation of these models.  Additionally, since the flutter speed predicted by the Wagner and Peters' models match, we can be reasonably confident that the Wagner unsteady aerodynamic model is also implemented correctly.

Time domain simulations may also be used in order to determine a system's stability.  To perform time domain simulations, an object representing the ordinary differential equations corresponding to the model may be generated using the [`get_ode`](@ref) function and then solved using the [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl) package.

```@example typical-section-stability
using DifferentialEquations

# models
aerodynamic_model = Peters{6}()
structural_model = TypicalSection()
coupled_model = couple_models(aerodynamic_model, structural_model)

# non-dimensional parameters
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
cd0 = 0
cm0 = 0
V = 1.0 # = U/(b*ωθ) (reduced velocity)

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
U = V*b*ωθ

# parameters
p_aero = [a, b, a0, α0, cd0, cm0]
p_stru = [kh, kθ, m, Sθ, Iθ]
p_additional = [U, ρ]
p = vcat(p_aero, p_stru, p_additional)

# initial states
u0_aero = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
u0_stru = [1.0, 0.0, 0.0, 0.0] # non-zero plunge degree of freedom
x0 = vcat(u0_aero, u0_stru)

# simulate from 0 to 10 seconds
tspan = (0.0, 100.0)

# construct ODE function
f = get_ode(coupled_model)

# construct ODE problem
prob = DifferentialEquations.ODEProblem(f, x0, tspan, p)

# solve ODE
sol = DifferentialEquations.solve(prob)

nothing #hide
```

We can then plot the solution using DifferentialEquations' built-in interface with the [Plots](https://github.com/JuliaPlots/Plots.jl) package.

```@example typical-section-stability
using Plots
pyplot()

plot(sol,
    vars = [7,8,9,10],
    xlabel = "t",
    ylabel = permutedims([
        "\$h\$",
        "\$\\theta\$",
        "\$\\dot{h}\$",
        "\$\\dot{\\theta}\$",
        ]),
    label = "",
    layout = (4, 1),
    size = (600,1200)
    )

savefig("typical-section-solution.svg") #hide

nothing #hide
```

![](typical-section-solution.svg)

We can also visualize and/or create animations of the solution geometry with the help of custom plot recipes.

```@example typical-section-stability
# Plot Recipes:
# - plot(model, dx, x, p, t)
# - plot(model, sol) # plots the solution geometry at index `sol.tslocation`
# - plot(model, sol, t) # plots the solution geometry at time `t`

# create animation
anim = @animate for t in range(tspan[1], tspan[2], length=200)
    plot(coupled_model, sol, t)
end

# save animation
gif(anim, "typical-section-simulation.gif")

nothing #hide
```

![](typical-section-simulation.gif)

Visualizing the solution geometry is especially helpful for identifying mode shapes.  For example, we can visualize the mode shape of the flutter mode for this system.

```@example typical-section-stability
# Plot Recipes:
# - plot(model, dx, x, p, t)
# - plot(model, sol) # plots the solution geometry at index `sol.tslocation`
# - plot(model, sol, t) # plots the solution geometry at time `t`

# create animation
anim = @animate for t in range(tspan[1], tspan[2], length=200)
    plot(coupled_model, sol, t)
end

# save animation
gif(anim, "typical-section-flutter-mode.gif")

nothing #hide
```

![](typical-section-flutter-mode.gif)

## Aeroelastic Analysis of the Goland Wing

![](goland-wing.png)

In this example, we demonstrate how to perform a three-dimensional aeroelastic analysis using the Goland wing, a low-aspect ratio prismatic metallic wing, which has been extensively used for validation.

```@setup goland-stability
using Plots
pyplot()
nothing #hide
```

The Goland wing is a cantilevered wing with a 20 ft span and 6 ft chord.  Its airfoil consists of a 4% thick parabolic arc.  There are two configurations for this wing, one with a tip store and one without.  The configuration we consider in this example is the configuration without a tip store.

The deflections of Goland wing are relatively small, so linear structural models are sufficient for modeling the wing's structure.  However, to demonstrate the capabilities of this package, we will use a nonlinear geometrically exact beam theory model.  

For the aerodynamics, we use a lifting line model which is capable of using a variety of 2D models to model section lift and moment coefficients.  While this type of aerodynamic model is inappropriate for this wing due to the wing's low aspect ratio, we use it so that we can obtain a better comparison between our results and the results of other aeroelastic analyses of the Goland wing performed with lifting line aerodynamics.

```@example goland-stability
using AerostructuralDynamics, GXBeam, DifferentialEquations, LinearAlgebra

# discretization
N = 8 # number of elements

# geometric properties
span = 6.096 # m (20 ft span)
chord = 1.8288 # m (6 ft chord)

# structural section properties
xea = 0.33*chord # m (elastic axis, from leading edge)
EIcc = 9.77e6 # N*m^2 (flat bending rigidity)
GJ = 0.99e6 # N*m^2 (torsional rigidity)
μ = 35.71 # kg/m (mass per unit length)
xcm = 0.43*chord # m (center of mass, from leading edge)
i11 = 8.64 # kg*m (moment of inertia about elastic axis)
i22 = 0.1*i11 # moment of inertia about beam y-axis
i33 = 0.9*i11 # moment of inertia about beam z-axis

# freestream properties
Vinf = vcat(1, 5:5:200) # m/s (velocity)
α = 0 # angle of attack

# aerodynamic section properties
xref = xea/chord # normalized reference location (relative to leading edge)
a = xref - 0.5 # normalized reference location (relative to semi-chord)
b = chord / 2 # m (semi-chord)
ρ = 1.02 # kg/m^3 (air density)
a0 = 0.85*(2*pi) # lift slope (for each section)
α0 = 0 # zero lift angle of attack (for each section)

# define geometry
xpt = range(0, 0, length=N+1) # point x-coordinates (in body frame)
ypt = range(0, span, length=N+1) # point y-coordinates (in body frame)
zpt = range(0, 0, length=N+1) # point z-coordinates (in body frame)
points = [[xpt[i],ypt[i],zpt[i]] for i = 1:N+1]
start = 1:N # starting point of each beam element
stop = 2:N+1 # ending point of each beam element
frames = fill([0 1 0; 1 0 0; 0 0 -1], N) # local to body frame transformation
compliance = fill(Diagonal([0, 0, 0, 1/GJ, 1/EIcc, 0]), N) # compliance matrix
xm2 = xea - xcm
mass = fill([
    μ 0 0 0 0 -μ*xm2;
    0 μ 0 0 0 0;
    0 0 μ μ*xm2 0 0;
    0 0 μ*xm2 i11 0 0;
    0 0 0 0 i22 0;
    -μ*xm2 0 0 0 0 i33], N) # mass matrix
assembly = GXBeam.Assembly(points, start, stop; frames, compliance, mass)

prescribed = Dict(
    # fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0,
        theta_z=0),
)

# prescribed point conditions
point_conditions = zeros(6, length(assembly.points))

# additional distributed loads
element_loads = zeros(6, length(assembly.elements))

# construct aerodynamic model
aerodynamic_model = LiftingLine{N}(Peters{6}())

# construct structural model
structural_model = GEBT(assembly, prescribed)

# define coupled model
model = couple_models(aerodynamic_model, structural_model)

# construct ODE function
f = get_ode(model)

# current time
t = 0.0

# eigenvalue/eigenvector storage
nev = 12*N
λ = zeros(ComplexF64, nev, length(Vinf))
Uλ = zeros(ComplexF64, nev, number_of_states(model), length(Vinf))
Vλ = zeros(ComplexF64, number_of_states(model), nev, length(Vinf))

# initial guess for state variables
x0 = zeros(number_of_states(model))

# loop through each velocity
for i = 1:length(Vinf)

    println(Vinf[i])

    # set parameters
    p_aero = get_parameters(aerodynamic_model; section_parameters =
        fill((a = a, b = b, a0 = a0, alpha0 = α0), N))

    p_stru = get_parameters(structural_model; assembly = assembly)

    p_additional = get_additional_parameters(model;
        rho = ρ,
        point_conditions = point_conditions,
        element_loads = element_loads,
        x = Vinf[i],
        v = 0,
        w = 0,
        p = 0,
        q = 0,
        r = 0,
        )

    p = vcat(p_aero, p_stru, p_additional)

    # find state variables corresponding to steady state operating conditions
    sol = solve(SteadyStateProblem(f, x0, p), SSRootfind())

    # use state variables from steady state operating conditions
    x = sol.x

    # calculate inputs corresponding to steady state operating conditions
    y = get_coupling_inputs(model, x, p, t)

    # perform linear stability analysis
    λi, Uλi, Vλi = get_eigen(model, x, y, p, t; nev)

    # correlate eigenvalues
    if i > 1
        # previous left eigenvector matrix
        Uλpi = Uλ[:,:,i-1]

        # current mass matrix
        Mi = get_mass_matrix(model, x, y, p, t)

        # use correlation matrix to correlate eigenmodes
        perm, corruption = AerostructuralDynamics.correlate_eigenmodes(Uλpi, Mi, Vλi)

        # re-arrange eigenmodes
        λi = λi[perm]
        Uλi = Uλi[perm,:]
        Vλi = Vλi[:,perm]
    end

    # save eigenvalues/eigenvectors
    λ[:,i] = λi
    Uλ[:,:,i] = Uλi
    Vλ[:,:,i] = Vλi

    # update initial guess for the state variables
    x0 .= x
end
```

To identify the flutter speed and frequency, we can plot the results.

```@example goland-stability
using Plots
pyplot()

sp1 = plot(
    xlim = (0, 200),
    xtick = 0:40:200,
    xlabel = "Velocity (m/s)",
    ylim = (0, 800),
    ytick = 0:100:800,
    ylabel = "Frequency (rad/s)",
    framestyle = :zerolines,
    titlefontsize = 14,
    guidefontsize = 14,
    legendfontsize = 11,
    tickfontsize = 11,
    legend = :topright,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true)

sp2 = plot(
    xlim = (0, 200),
    xtick = 0:40:200,
    xlabel = "Velocity (m/s)",
    ylim = (-80, 20),
    ytick = -80:20:20,
    ylabel = "Damping (1/s)",
    framestyle = :zerolines,
    titlefontsize = 14,
    guidefontsize = 14,
    legendfontsize = 11,
    tickfontsize = 11,
    legend = :topright,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=true)

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    if any(-80 .<= real.(λi) .<= 20)
        plot!(sp1, Vi, imag.(λi),
            label = "",
            color = i,
            markersize = 3,
            markerstrokewidth = 0,
            )
    end

end

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    if any(-80 .<= real.(λi) .<= 20)
        plot!(sp2, Vi,
            real.(λi),
            label = "",
            color = i,
            markersize = 3,
            markerstrokewidth = 0,
            )
    end
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800), show=true)
```

As predicted by this analysis, the flutter speed is 140 m/s and the flutter frequency is 69.7 rad/s.  These results are nearly identical to results by Palacios and Epureanu in "An Intrinsic Description of the Nonlinear Aeroelasticity of Very Flexible Wings" for a similar aeroelastic model.  We can therefore be reasonably confident that the models used for this example are implemented correctly.

```@example typical-section-stability
using DifferentialEquations

Vinf = 100.0

x0 = zeros(number_of_states(model))

# set parameters
p_aero = get_parameters(aerodynamic_model; section_parameters =
    fill((a = a, b = b, a0 = a0, alpha0 = α0), N))

p_stru = get_parameters(structural_model; assembly = assembly)

p_additional = get_additional_parameters(model;
    rho = ρ,
    point_conditions = point_conditions,
    element_loads = element_loads,
    x = 100.0,
    v = 0,
    w = 0,
    p = 0,
    q = 0,
    r = 0,
    )

p = vcat(p_aero, p_stru, p_additional)

# update initial state variables with steady state solution
x0 .= solve(SteadyStateProblem(f, x0, p), SSRootfind())

# simulate from 0 to 10 seconds
tspan = (0.0, 10.0)

# construct ODE function
f = get_ode(model)

# construct ODE problem
prob = DifferentialEquations.ODEProblem(f, x0, tspan, p)

# solve ODE
sol = DifferentialEquations.solve(prob)

nothing #hide
```

## Aeroelastic Analysis of a Highly Flexible Cantilever Wing

![](cantilever-wing.png)

In this example, we demonstrate how to perform a three-dimensional aeroelastic analysis using a highly flexible cantilever wing.  

```@setup cantilever-stability
using Plots
pyplot()
nothing #hide
```

The wing we are considering in this example was created by modifying Daedalus aircraft data and is therefore representative of a high-altitude long-endurance wing.  It has a 16 meter span (from root to tip) and a 1 meter chord.  To model the wing's aerodynamics, we will use a lifting line model.  To model the wing's structure, we will use a geometrically exact beam theory model.

```@example cantilever-stability
using AerostructuralDynamics, GXBeam, DifferentialEquations, LinearAlgebra

# discretization
N = 1 # number of elements

# geometric properties
span = 16 # m
chord = 1 # m

# structural section properties
GJ = 1e4 # N*m^2 (torsional rigidity)
EIcc = 2e4 # N*m^2 (flat bending rigidity)
EInn = 4e6 # N*m^2 (chord bending rigidity)
μ = 0.75 # kg/m (mass per unit length)
i11 = 0.1 # kg*m (moment of inertia about elastic axis)
i22 = 0.1*i11 # moment of inertia about beam y-axis
i33 = 0.9*i11 # moment of inertia about beam z-axis

# freestream properties
Vinf = vcat(0.1, 1:35) # m/s (velocity)
α = 2*pi/180 # angle of attack

# aerodynamic section properties
a = 0.0 # normalized reference location (relative to semi-chord)
b = chord / 2 # m (semi-chord)
ρ = 0.088 # kg/m^3 (air density)
a0 = 2*pi # lift slope (for each section)
α0 = 0 # zero lift angle of attack (for each section)

# define geometry
xpt = range(0, 0, length=N+1) # point x-coordinates (in body frame)
ypt = range(0, span, length=N+1) # point y-coordinates (in body frame)
zpt = range(0, 0, length=N+1) # point z-coordinates (in body frame)
points = [[xpt[i],ypt[i],zpt[i]] for i = 1:N+1]
start = 1:N # starting point of each beam element
stop = 2:N+1 # ending point of each beam element
frames = fill([0 1 0; 1 0 0; 0 0 -1], N) # local to body frame transformation
compliance = fill(Diagonal([0, 0, 0, 1/GJ, 1/EIcc, 1/EInn]), N) # compliance matrix
mass = fill(Diagonal([μ, μ, μ, i11, i22, i33]), N)
assembly = GXBeam.Assembly(points, start, stop; frames, compliance, mass)

prescribed = Dict(
    # fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0,
        theta_z=0),
)

# prescribed point conditions
point_conditions = zeros(6, length(assembly.points))

# additional distributed loads
element_loads = zeros(6, length(assembly.elements))

# construct aerodynamic model
aerodynamic_model = LiftingLine{N}(Peters{6}())

# construct structural model
structural_model = GEBT(assembly, prescribed)

# define coupled model
model = couple_models(aerodynamic_model, structural_model)

# construct ODE function
f = get_ode(model)

# current time
t = 0.0

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(model), length(Vinf))

# initial guess for state variables
x0 = zeros(number_of_states(model))

# loop through each velocity
for i = 1:length(Vinf)

    # previous left eigenvector matrix
    global Uλpi

    # set parameters
    p_aero = get_parameters(aerodynamic_model; section_parameters =
        fill((a = a, b = b, a0 = a0, alpha0 = α0), N))

    p_stru = get_parameters(structural_model; assembly = assembly)

    p_additional = get_additional_parameters(model;
        rho = ρ,
        point_conditions = point_conditions,
        element_loads = element_loads,
        x = Vinf[i]*cos(α),
        v = 0,
        w = -Vinf[i]*sin(α),
        p = 0,
        q = 0,
        r = 0,
        )

    p = vcat(p_aero, p_stru, p_additional)

    # find state variables corresponding to steady state operating conditions
    sol = solve(SteadyStateProblem(f, x0, p), SSRootfind())

    # use state variables from steady state operating conditions
    x = sol.x

    # calculate inputs corresponding to steady state operating conditions
    y = get_coupling_inputs(model, x, p, t)

    # perform linear stability analysis
    λi, Uλi, Vλi = get_eigen(model, x, y, p, t; nev = number_of_states(model))

    # correlate eigenvalues
    if i > 1
        # calculate mass matrix
        Mi = get_mass_matrix(model, x, y, p, t)

        # use correlation matrix to correlate eigenmodes
        perm, corruption = AerostructuralDynamics.correlate_eigenmodes(Uλpi, Mi, Vλi)

        # re-arrange eigenmodes
        λi = λi[perm]
        Uλi = Uλi[perm,:]
        Vλi = Vλi[:,perm]
    end

    # save eigenvalues
    λ[:,i] = λi

    # save previous left eigenvector matrix
    Uλpi = Uλi

    # update initial guess for the state variables
    x0 .= x
end

nothing #hide
```

We now plot the results predicted using each aerodynamic model.

```@example cantilever-stability
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
    xlim = (0, 35),
    xtick = 0:5:35,
    xlabel = "Velocity (m/s)",
    ylim = (0, 350),
    ytick = 0:50:350,
    ylabel = "Frequency (rad/s)",
    legend = :topright
    )

sp2 = plot(
    xlim = (0, 35),
    xtick = 0:5:35,
    xlabel = "Velocity (m/s)",
    ylim = (-12, 8),
    ytick = -12:4:8,
    ylabel = "Damping (1/s)",
    legend = :topleft
    )

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    if any(abs.(λi) .< 1e4)
        plot!(sp1, Vi, imag.(λi),
            label = "",
            color = i,
            markersize = 3,
            markerstrokewidth = 0,
            )
    end

end

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    if any(abs.(λi) .< 1e4)
        plot!(sp2, Vi,
            real.(λi)./abs.(λi)*100,
            label = "",
            color = i,
            markersize = 3,
            markerstrokewidth = 0,
            )
    end
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800), show=true)

nothing #hide
```

As predicted by this analysis, the flutter speed is 140 m/s and the flutter frequency is 69.7 rad/s.  These results are nearly identical to results by Palacios and Epureanu in "An Intrinsic Description of the Nonlinear Aeroelasticity of Very Flexible Wings" for a similar aeroelastic model.  We can therefore be reasonably confident that the models used for this example are implemented correctly.

## Aeroelastic Analysis of a Blended-Wing-Body Aircraft

In this example, we demonstrate how to perform a three-dimensional aeroelastic analysis using geometrically exact beam theory in combination with a rigid body dynamics model and various aerodynamic models.  

```@setup blended-wing-body-stability
using Plots
pyplot()
nothing #hide
```

```julia
using AerostructuralDynamics, GXBeam, NLsolve, LinearAlgebra

# discretization
NB = 6 # number of elements on body
NW = 12 # number of elements on wing
N = NB + NW # number of elements on blended wing body aircraft

# wing geometry
chord = [1.39, 0.55, 0.55] # m (chord lengths in -x-direction)
ybreak = [0.0, 0.89, 3.25] # m (span break locations in y-direction)
sweep = 30 * pi/180 # leading edge sweep

# stiffness properties
EA = [1.69e8, 1.55e8] # N (extension stiffness)
GJ = [2.25e6, 1.10e4] # N*m^2 (torsion stiffness)
EIcc = [7.5e5, 1.14e4] # N*m^2 (out-of-plane bending stiffness)
EInn = [3.5e7, 1.3e5] # N*m^2 (in-plane bending stiffness)

# inertial properties
μ = [50.0, 6.2] # kg/m (mass per unit length)
i11 = [4.5, 5.08e-3] # kg*m (rotational inertia per unit length)
i22 = [0.7, 5.0e-4] # kg*m (flat bending inertia per unit length)
i33 = [22, 4.63e-3] # kg*m (edge bending inertia per unit length)

# axis locations
xref = [0.6438, 0.4560, 0.4560] .* chord # reference axis location (from leading edge)
xcm = [0.6438, 0.4560, 0.4560] .* chord # center of gravity (from leading edge)

# body definition
xpt1 = range(0 + xref[1], ybreak[1]*tan(sweep) + xref[2]; length=NB+1)
ypt1 = range(0, ybreak[1]; length=NB+1)
zpt1 = range(0, 0; length=NB+1)
points1 = [[xpt1[i],ypt1[i],zpt1[i]] for i = 1:NB+1]
frames1 = fill( , NB)
compliance1 = fill(Diagonal([1/EA, 0, 0, 1/GJ, 1/EIcc, 1/EInn]), NB)
mass1 = fill(Diagonal([μ[1], μ[1], μ[1], i11[1], i22[1], i33[1]]), NB)
a1 = linterp(xref[1]/chord[1] - 0.5, xref[2]/chord[2] - 0.5, length=NB)
b1 = linterp(chord[1]/2, chord[2]/2, length=NB) # m (semi-chord)
a01 = fill(2*pi, NB) # lift slope (for each section)
α01 = fill(0.0, NB) # zero lift angle of attack (for each section)

# wing definition
xpt2 = range(ybreak[1]*tan(sweep) + xref[2], ybreak[2]*tan(sweep) + xref[3];
    length=NW+1)
ypt2 = range(ybreak[1], ybreak[2]; length=NW+1)
zpt2 = range(0, 0; length=NW+1)
points2 = [[xpt2[i],ypt2[i],zpt2[i]] for i = 1:NW+1]
frames2 = fill( , NW)
compliance2 = fill(Diagonal([1/EA[2], 0, 0, 1/GJ[2], 1/EIcc[2], 1/EInn[2]]), NW)
mass2 = fill(Diagonal([μ[2], μ[2], μ[2], i11[2], i22[2], i33[2]]), NW)
a2 = linterp(xref[2]/chord[2] - 0.5, xref[3]/chord[3] - 0.5; length = NW)
b2 = linterp(chord[2]/2, chord[3]/2; length = NW) # m (semi-chord)
a02 = fill(2*pi, NW) # lift slope (for each section)
α02 = fill(0.0, NW) # zero lift angle of attack (for each section)

# define beam assembly
points = vcat(points1, points2[2:end])
frames = vcat(frames1, frames2)
start = 1:N # starting point of each beam element
stop = 2:N+1 # ending point of each beam element
assembly = GXBeam.Assembly(points, start, stop; frames, compliance, mass)

# define aerodynamic section properties
a = vcat(a1, a2)
b = vcat(b1, b2)
a0 = vcat(a01, a02)
α0 = vcat(α01, α02)
ρ = 1.02 # kg/m^3 (air density)

# boundary condition initialization
prescribed = Dict(
    # fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0,
        theta_z=0),
)

# construct aerodynamic model
aerodynamic_model = LiftingLine{N}(Peters{6}())

# construct structural model
structural_model = GEBT(assembly, prescribed)

# rigid body dynamics model
dynamics_model = RigidBody()

# define simulation models
model = couple_models(aerodynamic_model, structural_model, dynamics_model)

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(model), length(Vinf))

x0 = zeros(number_of_states(model))

# loop through each velocity
for i = 1:length(Vinf)

    println("Vinf: ", Vinf[i])

    # set state variables, parameters, and current time
    p_aero = vcat([(a[i], b[i], a0[i], α0[i]) for i = 1:N]...)
    p_stru = set_parameters(structural_model, assembly)
    p_additional = vcat(-Vinf[i]*cos(α), 0, -Vinf[i]*sin(α), ρ,
        set_inputs(structural_model, assembly; prescribed=prescribed))
    p = vcat(p_aero, p_stru, p_additional)
    t = 0

    # find state variables corresponding to steady state operating conditions
    fresid = x -> get_rates(model, x, get_coupling_inputs(model, x, p, t), p, t)
    sol = nlsolve(fresid, x0)
    x = sol.zero

    # calculate the inputs corresponding to steady state operating conditions
    y = get_coupling_inputs(model, x, p, t)

    # calculate the mass matrix corresponding to steady state operating conditions
    M = get_mass_matrix(model, x, y, p, t)

    # calculate the jacobian corresponding to steady state operating conditions
    J = get_state_jacobian(model, x, y, p, t)

    # solve the generalized eigenvalue problem
    λ[:,i] = sort(eigvals(J, M), by=LinearAlgebra.eigsortby)

    # update initial guess for the state variables
    x0 .= x
end
```

We now plot the results predicted using each aerodynamic model.

```julia
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
    xlim = (0, 200),
    xtick = 0:40:200,
    xlabel = "Velocity (m/s)",
    ylim = (0, 1000),
    ytick = 0:100:1000,
    ylabel = "Frequency (rad/s)",
    legend = :topright
    )

sp2 = plot(
    xlim = (0, 200),
    xtick = 0:40:200,
    xlabel = "Velocity (m/s)",
    ylim = (-80, 20),
    ytick = -80:20:20,
    ylabel = "Damping (1/s)",
    legend = :topleft
    )

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    scatter!(sp1, Vi, imag.(λi),
        label = "",
        color = 1,
        markersize = 3,
        markerstrokewidth = 0,
        )
end

for i = 1:size(λ, 1)

    Vi = Vinf[:]
    λi = λ[i,:]

    scatter!(sp2, Vi,
        real.(λi),
        label = "",
        color = 1,
        markersize = 3,
        markerstrokewidth = 0,
        )
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800), show=true)

```
