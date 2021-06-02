# Examples

## Aeroelastic Analysis of the Typical Section Model

In this example, we demonstrate how to perform a two-dimensional aeroelastic analysis using a typical section model with two degrees of freedom.

![](typical-section.svg)

The equations of motion for this model are
```math
m \left(\ddot{h}+b x_\theta \ddot{\theta} \right) + k_h h = -L \\
I_P \ddot{\theta} + m b x_\theta \ddot{h} + k_\theta = M_{\frac{1}{4}} + b \left( \frac{1}{2} + a \right) L
```
where ``a`` is the normalized distance from the semichord to the reference point, ``b`` is the semichord length, ``k_h`` is the linear spring constant, ``k_\theta`` is the torsional spring constant, ``m`` is the mass per unit span, ``x_\theta`` is the distance to the center of mass from the reference point, ``I_P`` is the moment of inertia about the reference point, ``L`` is the lift per unit span, and ``M_\frac{1}{4}`` is the quarter-chord moment per unit span.

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

```@example typical-section-stability
using AerostructuralDynamics, LinearAlgebra

# reduced velocity range
V = range(0, 3.1, length=5000) # (reduced velocity)

# non-dimensional parameters
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Ip/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a

# chosen dimensional parameters
b = 1
ρ = 1
ωθ = 1

# derived dimensional parameters
m = μ*ρ*pi*b^2
Ip = r2*m*b^2
ωh = σ*ωθ
kh = m*ωh^2
kθ = Ip*ωθ^2

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
        p_aero = [a, b, U[i], ρ]
        p_stru = [a, b, kh, kθ, m, xθ, Ip]
        p = vcat(p_aero, p_stru)

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

```@example typical-section-stability
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

savefig(p1, "typical-section-stability.svg") #hide

nothing #hide
```

![](typical-section-stability.svg)

The same analysis and results are presented by Hodges and Pierce in "Introduction to Structural Dynamics and Aeroelasticity" for the steady state and Peters' Finite State models.  The results shown here match with those provided by Hodges and Pierce, thus validating our implementation of these models.
