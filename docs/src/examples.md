# Examples

## Two-Dimensional Aeroelastic Analysis

In this example, we demonstrate how to perform a two-dimensional aeroelastic analysis using a typical section model with two degrees of freedom.

![](typical-section.svg)

The equations of motion for this model are
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M
```
where ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_θ`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M`` is the moment per unit span about the reference point.

We use the non-dimensional parameters
```math
a = -1/5 \quad e = -1/10  \\
r^2 = \frac{I_P}{m b^2} \quad \sigma = \frac{\omega_h}{\omega_\theta} \\
\mu = \frac{m}{\rho_\infty \pi b^2} \quad V = \frac{U}{b \omega_\theta}
```
where ``a`` is the normalized distance from the semichord to the reference point, ``e`` is the normalized distance from the semichord to the center of mass, and ``\omega_h`` and ``\omega_\theta`` are the uncoupled natural frequencies.
```math
\omega_h = \sqrt{\frac{k_h}{m}} \quad \omega_\theta = \sqrt{\frac{k_\theta}{I_P}}
```

We perform aeroelastic analyses using a variety of aerodynamic models in order to compare the various models.

```@setup two-dimensional-stability
using Plots
pyplot()
nothing #hide
```

```@example two-dimensional-stability
using AerostructuralDynamics, LinearAlgebra

# reduced velocity range
V = range(0, 3.1, length=5000) # (reduced velocity)

# non-dimensional parameters
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle

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

# eigenvalue storage
λ = Vector{Matrix{ComplexF64}}(undef, length(aerodynamic_models))

# loop through each aerodynamic model
for (ia, aerodynamic_model) in enumerate(aerodynamic_models)

    # combined models
    models = (aerodynamic_model, structural_model)

    # eigenvalue storage
    λ[ia] = zeros(ComplexF64, number_of_states(models), length(V))

    # loop through each reduced frequency
    for i = 1:length(V)
        # state variables
        u_aero = zeros(number_of_states(aerodynamic_model))
        u_stru = zeros(number_of_states(structural_model))
        u = vcat(u_aero, u_stru)

        # parameters
        p_aero = [a, b, ρ, a0, α0]
        p_stru = [kh, kθ, m, Sθ, Iθ]
        p_input = U[i]
        p = vcat(p_aero, p_stru, p_input)

        # time
        t = 0.0

        # calculate inputs
        y = get_inputs(models, u, p, t)

        # mass matrix
        M = get_mass_matrix(models, u, y, p, t)

        # jacobian
        J = get_state_jacobian(models, u, y, p, t)

        # solve generalized eigenvalue problem
        λ[ia][:,i] = sort(eigvals(J, M), by=LinearAlgebra.eigsortby)
    end
end

nothing #hide
```

We now plot the results predicted using each aerodynamic model.

```@example two-dimensional-stability
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
    title = "Non-Dimensional Frequency",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    legend = :topright
    )

sp2 = plot(
    title = "Non-Dimensional Damping",
    xlim = (0,3.1),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (-0.7, 0.605),
    ytick = -0.6:0.2:0.6,
    ylabel = "\$ \\frac{Γ}{\\omega_\\theta} \$",
    legend = :topleft
    )

labels = ["Steady", "Quasi-Steady", "Wagner", "Peters (N=6)"]

for ia = 1:length(aerodynamic_models)

    scatter!(sp1, V, imag.(λ[ia][1,:])/ωθ,
        label = labels[ia],
        color = ia,
        markersize = 1,
        markerstrokewidth = 0,
        )

    for i = 2:size(λ[ia], 1)
        scatter!(sp1, V, imag.(λ[ia][i,:])/ωθ,
            label = "",
            color = ia,
            markersize = 1,
            markerstrokewidth = 0,
            )
    end

    scatter!(sp2, V, real.(λ[ia][1,:])/ωθ,
        label = labels[ia],
        color = ia,
        markersize = 1,
        markerstrokewidth = 0,
        )

    for i = 2:size(λ[ia], 1)
        scatter!(sp2, V, real.(λ[ia][i,:])/ωθ,
            label = "",
            color = ia,
            markersize = 1,
            markerstrokewidth = 0,
            )
    end
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800))

savefig(p1, "two-dimensional-stability.svg") #hide

nothing #hide
```

![](two-dimensional-stability.svg)

The same analysis and results are presented by Hodges and Pierce in "Introduction to Structural Dynamics and Aeroelasticity" for the steady state and Peters' Finite State aerodynamic models.  The results shown here match with those provided by Hodges and Pierce, thus validating our implementation of these models.

## Three Dimensional Aeroelastic Analysis

In this example, we demonstrate how to perform a three-dimensional aeroelastic analysis using geometrically exact beam theory in combination with various aerodynamic models.

```@setup three-dimensional-stability
using Plots
pyplot()
nothing #hide
```

```@example three-dimensional-stability
using AerostructuralDynamics, GXBeam, NLsolve, LinearAlgebra

# aerodynamic properties
Vinf = range(1, 30, length=200) # m/s (velocity)
α = 2*pi/180 # angle of attack
ρ = 0.088 # kg/m^3 (air density)
a0 = 2*pi # lift slope (for each section)
α0 = 0 # zero lift angle of attack (for each section)

# structural properties
a = 0 # normalized reference axis location relative to the semichord
b = 16 # m (span)
c = 1 # m (chord)
GJ = 1e4 # N*m^2 (torsional rigidity)
EIcc = 4e6 # N*m^2 (chord bending rigidity)
EInn = 2e4 # N*m^2 (flat bending rigidity)
μ = 0.75 # kg/m (mass per unit span)
i11 = 0.1 # kg*m (rotational inertia per unit span)

# discretization
N = 8 # number of elements

# geometry initialization
x = range(0, b, length=N+1) # point x-coordinates
y = range(0, 0, length=N+1) # point y-coordinates
z = range(0, 0, length=N+1) # point z-coordinates
points = [[x[i],y[i],z[i]] for i = 1:N+1]
start = 1:N # starting point of each beam element
stop = 2:N+1 # ending point of each beam element
compliance = fill(Diagonal([0, 0, 0, 1/GJ, 1/EIcc, 1/EInn]), N) # compliance matrix
mass = fill(Diagonal([μ, μ, μ, i11, 0, 0]), N) # mass matrix
assembly = GXBeam.Assembly(points, start, stop; compliance = compliance,
    mass = mass)

# boundary condition initialization
prescribed = Dict(
    # fixed left edge
    1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0,
        theta_z=0),
)

# distributed load initialization
distributed = Dict()
for i = 1:N
    # distributed load on each beam element
    distributed[i] = GXBeam.DistributedLoads(assembly, i)
end

# structural system initialization
ipoint = [1] # points at which prescribed conditions are applied
static = false # static simulation?
system = GXBeam.System(assembly, ipoint, static)

# construct aerodynamic model
aerodynamic_model = LiftingLine{N}(Wagner())

# construct structural model
structural_model = GEBT(system, assembly, prescribed, distributed)

# define simulation models
models = (aerodynamic_model, structural_model)

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(models), length(Vinf))

# state variable initial guess
u0 = zeros(number_of_states(models))

# loop through each velocity
for i = 1:length(Vinf)

    println("Vinf: ", Vinf[i])

    # set parameters
    p_aero = vcat(fill([a, c/2, ρ, a0, α0], N)...)
    p_stru = Float64[]
    p_coupled = [0, -Vinf[i]*cos(α), Vinf[i]*sin(α)] # freestream velocity in global frame
    p = vcat(p_aero, p_stru, p_coupled)

    # set time
    t = 0.0

    # find state variables corresponding to steady state operating conditions
    fresid = u -> get_rates(models, u, get_inputs(models, u, p, t), p, t)
    sol = nlsolve(fresid, u0)

    # set states to their steady-state values
    u = sol.zero

    # calculate the inputs corresponding to steady state operating conditions
    y = get_inputs(models, u, p, t)

    # calculate the mass matrix corresponding to steady state operating conditions
    M = get_mass_matrix(models, u, y, p, t)

    # calculate the jacobian corresponding to steady state operating conditions
    J = get_state_jacobian(models, u, y, p, t)

    # solve the generalized eigenvalue problem
    λ[:,i] = sort(eigvals(J, M), by=LinearAlgebra.eigsortby)

    # update initial guess for the state variables
    u0 .= u
end

nothing #hide
```

We now plot the results predicted using each aerodynamic model.

```@example three-dimensional-stability
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
    xlim = (1, 30),
    xtick = vcat(1, 5:5:30),
    xlabel = "Velocity (m/s)",
    ylim = (0, 350),
    ytick = 0.0:50:350,
    ylabel = "Frequency (rad/s)",
    legend = :topright
    )

sp2 = plot(
    xlim = (1, 30),
    xtick = vcat(1, 5:5:30),
    xlabel = "Velocity (m/s)",
    ylim = (-12, 8),
    #ytick = -12:4:8,
    ylabel = "Damping Ratio %",
    legend = :topleft
    )

for i = 1:size(λ, 1)
    scatter!(sp1, Vinf, imag.(λ[i,:]),
        label = "",
        color = 1,
        markersize = 3,
        markerstrokewidth = 0,
        )
end

for i = 1:size(λ, 1)
    scatter!(sp2, Vinf,
        #real.(λ[i,:]),
        real.(λ[i,:])./abs.(λ[i,:])*100,
        label = "",
        color = 1,
        markersize = 3,
        markerstrokewidth = 0,
        )
end

p1 = plot(sp1, sp2, layout = (2, 1), size = (600, 800))

# savefig(p1, "three-dimensional-stability.svg") #hide

nothing #hide
```
