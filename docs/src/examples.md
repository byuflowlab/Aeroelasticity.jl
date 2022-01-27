# Examples

## Aeroelastic Analysis of a Typical Section

```@setup typical-section
using Plots
pyplot()
nothing #hide
```

In this example, we perform a two-dimensional aeroelastic analysis for a two degree of freedom typical section model.

![](typical-section.svg)

The equations of motion for this model are
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_\theta \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M
```
where ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_θ`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M`` is the moment per unit span about the reference point.  To learn more about how this model is implemented in `AerostructuralDynamics`, see [`TypicalSection`](@ref).

AerostructuralDynamics has a number of pre-implemented aerodynamic models which may be used to model the aerodynamics of the typical section model.  These models include a steady-state thin airfoil theory model (see [`Steady`](@ref)), a quasi-steady thin airfoil theory model (see [`QuasiSteady`](@ref)), an unsteady aerodynamic model based on Wagner's function (see [`Wagner`](@ref)), and Peters' finite state aerodynamic model (see [`Peters`](@ref)).  We will perform aeroelastic analyses using each of these models and compare the results.

The non-dimensional parameters we will use are
```math
a = -1/5 \quad e = -1/10  \\
r^2 = \frac{I_P}{m b^2} \quad \sigma = \frac{\omega_h}{\omega_\theta} \\
\mu = \frac{m}{\rho_\infty \pi b^2} \quad V = \frac{U}{b \omega_\theta}
```
where ``a`` is the normalized distance from the semichord to the reference point, ``e`` is the normalized distance from the semichord to the center of mass, and ``\omega_h`` and ``\omega_\theta`` are the uncoupled natural frequencies, defined as
```math
\omega_h = \sqrt{\frac{k_h}{m}} \quad \omega_\theta = \sqrt{\frac{k_\theta}{I_P}}
```

```@example typical-section-stability
using AerostructuralDynamics, DifferentialEquations, LinearAlgebra

# define non-dimensional parameters
V = range(1e-6, 3.1, length=1000) # = U/(b*ωθ) (reduced velocity)
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a # distance from center of mass to reference point
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
cd0 = 0 # drag coefficient
cm0 = 0 # moment coefficient

# choose dimensional parameters
b = 1 # semichord
ρ = 1 # air density
ωθ = 1 # pitch natural frequency
c = 343 # air speed of sound

# calculate dimensionalized parameters
U = V*b*ωθ # freestrean velocity
m = μ*ρ*pi*b^2 # mass
Sθ = m*xθ*b # mass imbalance
Iθ = r2*m*b^2 # inertia
ωh = σ*ωθ # plunge natural frequency
kh = m*ωh^2 # plunge spring constant
kθ = Iθ*ωθ^2 # pitch spring constant

# define aerodynamic models
aerodynamic_models = (Steady(), QuasiSteady(), Wagner(), Peters())

# initialize eigenvalue/eigenvector storage
λ = Vector{Matrix{ComplexF64}}(undef, length(models))
Uλ = Vector{Array{ComplexF64,3}}(undef, length(models))
Vλ = Vector{Array{ComplexF64,3}}(undef, length(models))

# perform an analysis for each aerodynamic model
for imodel = 1:length(aerodynamic_models)

    # define coupled model
    model = assemble_model(;
        aerodynamic_model = aerodynamic_models[imodel],
        structural_model = Section())

    # define initial guess for equilibrium states
    x0 = assemble_states(model)

    # eigenvalue/eigenvector storage
    nλ = number_of_states(models[imodel])
    λ[imodel] = zeros(ComplexF64, nλ, length(V))
    Uλ[imodel] = zeros(ComplexF64, nλ, nλ, length(V))
    Vλ[imodel] = zeros(ComplexF64, nλ, nλ, length(V))

    # loop through each reduced frequency
    for i = 1:length(V)

        # define aerodynamic parameters
        aerodynamic_parameters = (; a = a, b = b, a0 = a0, alpha0 = α0, cd0 = cd0, cm0 = cm0)

        # define structural parameters
        structural_parameters = (; kh = kh, ktheta = kθ, m = m, Stheta = Sθ, Itheta = Iθ)

        # define additional parameters
        additional_parameters = (; U = U[i], rho = ρ, c = c)

        # define parameter vector
        p = assemble_parameters(model;
            aerodynamic_parameters = aerodynamic_parameters,
            structural_parameters = structural_parameters,
            additional_parameters = additional_parameters)

        # find equilibrium point
        x = solve(SteadyStateProblem(f, x0, p))

        # linearize about equilibrium point
        K, M = linearize(model, x, p)

        # perform linear stability analysis
        λi, Uλi, Vλi = get_eigen(model, K, M)

        # correlate eigenvalues
        if i > 1
            # previous left eigenvector matrix
            Uλpi = Uλ[imodel][:,:,i-1]

            # use correlation matrix to correlate eigenmodes
            perm, corruption = correlate_eigenmodes(Uλpi, M, Vλi)

            # re-arrange eigenmodes
            λi = λi[perm]
            Uλi = Uλi[perm,:]
            Vλi = Vλi[:,perm]
        end

        # save eigenvalues/eigenvectors
        λ[imodel][:,i] = λi
        Uλ[imodel][:,:,i] = Uλi
        Vλ[imodel][:,:,i] = Vλi
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
    titlefontsize = 10,
    guidefontsize = 10,
    legendfontsize = 8,
    tickfontsize = 9,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid=false
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
    titlefontsize = 10,
    guidefontsize = 10,
    legendfontsize = 8,
    tickfontsize = 9,
    legend = :topleft,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    minorgrid = false
    )

labels = ["Steady", "Quasi-Steady", "Wagner", "Peters (N=6)"]

for ia = 1:length(models)

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

p1 = plot(sp1, sp2, layout = (1, 2), size = (800, 300))

savefig(p1, "typical-section-stability.svg") #hide

nothing #hide
```

![](typical-section-stability.svg)

Using the `Wagner` or `Peters` aerodynamic models yields a flutter reduced velocity around 2.2, while the `Steady` and `QuasiSteady` aerodynamic models predict significantly lower flutter velocities.  The aerodynamic state variables of the `Wagner` and `Peters` models allows these models to capture the impact of vortex shedding on the lift and drag of the profile, therefore we can expect these models to yield more accurate results than the `Steady` and `QuasiSteady` models.

The non-dimensional parameters we use for this example match those used by Hodges and Pierce in "Introduction to Structural Dynamics and Aeroelasticity".  Hodges and Pierce performed the analysis using a steady-state model and Peter's finite state model with six state variables.   The results presented here for the steady-state and Peters' finite state models match the results presented by Hodges and Pierce in "Introduction to Structural Dynamics and Aeroelasticity", which validates our implementation of these models.  Additionally, since the flutter speed predicted by the `Wagner` and `Peters` models match, we can be reasonably confident that the Wagner unsteady aerodynamic model is also implemented correctly.

## Time Domain Simulation of a Typical Section

Time domain simulations may also be used in order to determine a system's stability.  To perform time domain simulations, an object representing the ordinary differential equations corresponding to the model may be generated using the [`get_ode`](@ref) function and then solved using the [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl) package.  For this example we use the same parameters as in the previous example.

```@example typical-section-simulation
using AerostructuralDynamics, DifferentialEquations, LinearAlgebra

# define non-dimensional parameters
V = range(1e-6, 3.1, length=1000) # = U/(b*ωθ) (reduced velocity)
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a # distance from center of mass to reference point
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
cd0 = 0 # drag coefficient
cm0 = 0 # moment coefficient

# choose dimensional parameters
b = 1 # semichord
ρ = 1 # air density
ωθ = 1 # pitch natural frequency
c = 343 # air speed of sound

# calculate dimensionalized parameters
U = V*b*ωθ # freestrean velocity
m = μ*ρ*pi*b^2 # mass
Sθ = m*xθ*b # mass imbalance
Iθ = r2*m*b^2 # inertia
ωh = σ*ωθ # plunge natural frequency
kh = m*ωh^2 # plunge spring constant
kθ = Iθ*ωθ^2 # pitch spring constant

# model
model = PetersSection(6)

# reduced velocity
V = 1.0 # = U/(b*ωθ)

# dimensionalized velocity
U = V*b*ωθ

# parameters
p_aero = [a, b, a0, α0, cd0, cm0]
p_stru = [kh, kθ, m, Sθ, Iθ]
p_additional = [U, ρ, c]
p = vcat(p_aero, p_stru, p_additional)

# initial states
u0_aero = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
u0_stru = [0.5, 0.0, 0.0, 0.0] # non-zero plunge degree of freedom
x0 = vcat(u0_aero, u0_stru)

# simulate from 0 to 100 seconds
tspan = (0.0, 100.0)

# construct ODE function
f = get_ode(model)

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

For aeroelastic models based on a typical section, we can also easily visualize the section's behavior.

```@example typical-section-stability

# animation parameters
a = -1/5
b = 0.5

# create animation
anim = @animate for t in range(tspan[1], tspan[2], length=200)
    h, θ = sol(t, idxs=7:8)
    xplot, yplot = section_coordinates(h, θ; a, b)
    plot(xplot, yplot;
        framestyle = :origin,
        grid = :false,
        xlims = (-1.0, 1.0),
        ylims = (-0.75, 0.75),
        aspect_ratio = 1.0,
        label = "t = $(round(t, digits=1))")
end

gif(anim, "typical-section-simulation.gif")

nothing #hide
```

![](typical-section-simulation.gif)

## Aeroelastic Analysis of the Goland Wing

![](goland-wing.png)

In this example, we demonstrate how to perform a three-dimensional aeroelastic analysis using the Goland wing, a low-aspect ratio prismatic metallic wing, which has been extensively used for validation.

```@setup goland-stability
using Plots
pyplot()
nothing #hide
```

The Goland wing is a cantilevered wing with a 20 ft span and 6 ft chord. Its airfoil consists of a 4% thick parabolic arc. There are two configurations for this wing, one with a tip store and one without. The configuration we consider in this example is the configuration without a tip store.

The deflections of Goland wing are relatively small, so they may be accurately modeled with a linear or nonlinear structural model.  For the purposes of this example we will use a low-order geometrically exact beam theory model with eight beam elements.  

To model the aerodynamics of the wing, we use a lifting line model.  While this type of aerodynamic model is inappropriate for this wing due to the wing's low aspect ratio, we use it so that we can obtain a better comparison between our results and the results of other aeroelastic analyses of the Goland wing performed with lifting line aerodynamics.

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

# Goland Wing Geometry
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