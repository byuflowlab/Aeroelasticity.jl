# Examples

## Typical Section Model with Steady Aerodynamics

```@example typical-section-steady
using AerostructuralDynamics, LinearAlgebra

# dimensionless parameters
a = -1/5 # center of mass normalized location
e = -1/10 # reference point normalized location
μ = 20 # = m/(ρ*pi*b^2)
r2 = 6/25 # = Ip/(m*b^2)
σ = 2/5 # = ωh/ωθ
V = range(1e-6, 3, length=100) # = U/(b*ωθ)

# chosen dimensional parameters
b = 1
ρ = 1
ωθ = 1

# dimensionalized parameters
m = μ*ρ*pi*b^2
Ip = r2*m*b^2
ωh = σ*ωθ
kh = m*ωh^2
kθ = Ip*ωθ^2
xθ = abs(e - a)
U = V*b*ωθ

# structural model
models = AerostructuralDynamics.CoupledTypicalSection()

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(models), length(V))

for i = 1:length(V)

    # state variables
    u = zeros(number_of_states(models))

    # inputs
    y = []

    # parameters
    p = [U[i], ρ, a, b, kh, kθ, m, xθ, Ip]

    # time
    t = 0.0

    # mass matrix
    M = get_mass_matrix(models, u, y, p, t)

    # jacobian
    J = get_state_jacobian(models, u, y, p, t)

    # solve generalized eigenvalue problem
    λ[:,i] = eigvals(J, M)

    # scale eigenvalues
    λ[:,i] .*= b/U[i]

end

nothing #hide
```

```@example typical-section-flutter
using Plots
pyplot()

p1 = scatter(V, imag.(λ')/(b*ωθ),
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    )

p2 = scatter(V, real.(λ')/(ωθ),
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (-0.7, 0.605),
    ytick = -0.6:0.2:0.6,
    ylabel = "\$ \\frac{Γ}{\\omega_\\theta} \$",
    )

plot(p1, p2, layout = (1, 2))

nothing #hide
```

```@example typical-section-steady
using AerostructuralDynamics, LinearAlgebra

# dimensionless parameters
a = -1/5 # center of mass normalized location
e = -1/10 # reference point normalized location
μ = 20 # = m/(ρ*pi*b^2)
r2 = 6/25 # = Ip/(m*b^2)
σ = 2/5 # = ωh/ωθ
V = range(1e-6, 3, length=100)

# dimensionalized parameters
b = 0.5
ρ = 1
m = ρ*pi*b^2
ωθ = 1
ωh = ωθ*σ
kh = m*ωh^2
kθ = m*ωθ^2
xθ = e - a
Ip = r2*m*b^2
U = V*b*ωθ

# aerodynamic model
aero = Steady()

# structural model
stru = TypicalSection()

# combined models
models = (aero, stru)

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(models), length(V))

for i = 1:length(V)

    # state variables
    u_aero = zeros(number_of_states(aero))
    u_stru = zeros(number_of_states(stru))
    u = vcat(u_aero, u_stru)

    # parameters
    p_aero = [b, U[i], ρ]
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
    λ[:,i] = eigvals(J, M)

end

nothing #hide
```

```@example typical-section-flutter
using Plots
pyplot()

p1 = scatter(V, imag.(λ')/(b*ωθ),
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    )

p2 = scatter(V, real.(λ')/(ωθ),
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (-0.7, 0.605),
    ytick = -0.6:0.2:0.6,
    ylabel = "\$ \\frac{Γ}{\\omega_\\theta} \$",
    )

plot(p1, p2, layout = (1, 2))

nothing #hide
```


## Typical Section Flutter Analysis

```@example typical-section-flutter
using AerostructuralDynamics, LinearAlgebra

# number of aerodynamic states
N = 6

# dimensionless parameters
a = -1/5 # center of mass normalized location
e = -1/10 # reference point normalized location
μ = 20 # = m/(ρ*pi*b^2)
r2 = 6/25 # = Ip/(m*b^2)
σ = 2/5 # = ωh/ωθ
V = range(1e-6, 3, length=100)

# dimensionalized parameters
b = 0.5
ρ = 1
m = ρ*pi*b^2
ωθ = 1
ωh = ωθ*σ
kh = m*ωh^2
kθ = m*ωθ^2
xθ = e - a
Ip = r2*m*b^2
U = V*b*ωθ

# aerodynamic model
aero = PetersFiniteState{N}()

# structural model
stru = TypicalSection()

# combined models
models = (aero, stru)

# eigenvalue storage
λ = zeros(ComplexF64, number_of_states(models), length(V))

for i = 1:length(V)

    # state variables
    u_aero = zeros(N)
    u_stru = zeros(4)
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
    λ[:,i] = eigvals(J, M)

end

nothing #hide
```

```@example typical-section-flutter
using Plots
pyplot()

p1 = scatter(V, imag.(λ')/(b*ωθ),
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    )

p2 = scatter(V, real.(λ')/(ωθ),
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (-0.7, 0.605),
    ytick = -0.6:0.2:0.6,
    ylabel = "\$ \\frac{Γ}{\\omega_\\theta} \$",
    )

plot(p1, p2, layout = (1, 2))

nothing #hide
```
