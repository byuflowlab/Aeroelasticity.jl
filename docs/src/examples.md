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
V = range(0, 3, length=1000)

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
    λ[:,i] = sort(eigvals(J, M), by=LinearAlgebra.eigsortby)

end

nothing #hide
```

```@example typical-section-steady
using Plots
pyplot()

p1 = scatter(V, imag.(λ')/ωθ, label="",
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    )

p2 = scatter(V, real.(λ')/ωθ, label="",
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

## Typical Section Model with Quasi-Steady Aerodynamics

```@example typical-section-quasisteady
using AerostructuralDynamics, LinearAlgebra

# dimensionless parameters
a = -1/5 # center of mass normalized location
e = -1/10 # reference point normalized location
μ = 20 # = m/(ρ*pi*b^2)
r2 = 6/25 # = Ip/(m*b^2)
σ = 2/5 # = ωh/ωθ
V = range(0, 3, length=1000)

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

# aerodynamic model
aero = QuasiSteady()

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
    λ[:,i] = sort(eigvals(J, M), by=LinearAlgebra.eigsortby)

end

nothing #hide
```

```@example typical-section-quasisteady
using Plots
pyplot()

p1 = scatter(V, imag.(λ')/ωθ, label="",
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    )

p2 = scatter(V, real.(λ')/ωθ, label="",a
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

## Typical Section Model with Peter's Finite State Aerodynamics

```@example typical-section-peters
using AerostructuralDynamics, LinearAlgebra

# dimensionless parameters
a = -1/5 # center of mass normalized location
e = -1/10 # reference point normalized location
μ = 20 # = m/(ρ*pi*b^2)
r2 = 6/25 # = Ip/(m*b^2)
σ = 2/5 # = ωh/ωθ
V = range(0, 3, length=1000)

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

# aerodynamic model
aero = PetersFiniteState{6}()

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
    λ[:,i] = sort(eigvals(J, M), by=LinearAlgebra.eigsortby)

end

nothing #hide
```

```@example typical-section-peters
using Plots
pyplot()

p1 = scatter(V, imag.(λ')/ωθ, label="",
    xlim = (0,3),
    xtick = 0.0:0.5:3.0,
    xlabel = "\$ \\frac{U}{b \\omega_\\theta} \$",
    ylim = (0, 1.05),
    ytick = 0.0:0.2:1.0,
    ylabel = "\$ \\frac{\\Omega}{\\omega_\\theta} \$",
    )

p2 = scatter(V, real.(λ')/ωθ, label="",a
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
