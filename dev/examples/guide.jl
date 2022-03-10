using Aeroelasticity

model = assemble_model(
    aerodynamic_model = Peters{4}(),
    structural_model = Section())

# non-dimensional parameters
V = 1.0 # = U/(b*ωθ) (reduced velocity)
a = -1/5 # reference point normalized location
e = -1/10 # center of mass normalized location
μ = 20 # = m/(ρ*pi*b^2) (mass ratio)
r2 = 6/25 # = Iθ/(m*b^2) (radius of gyration about P)
σ = 2/5 # = ωh/ωθ (natural frequency ratio)
xθ = e - a # normalized distance from the reference point to the center of mass
a0 = 2*pi # lift curve slope
α0 = 0 # zero lift angle
cd0 = 0 # zero lift drag coefficient
cm0 = 0 # zero lift moment coefficient

# chosen dimensional parameters
b = 0.5 # semi-chord
ρ = 1 # air density
ωθ = 1 # pitch natural frequency
c = 343 # air speed of sound

# dimensionalized parameters
U = V*b*ωθ # velocity
m = μ*ρ*pi*b^2 # mass
Sθ = m*xθ*b # structural imbalance
Iθ = r2*m*b^2 # moment of inertia
ωh = σ*ωθ # plunge natural frequency
kh = m*ωh^2 # linear spring constant
kθ = Iθ*ωθ^2 # torsional spring constant

p = assemble_parameters(model;
    aerodynamic_parameters = (a=a, b=b, a0=a0, alpha0=α0, cd0=cd0, cm0=cm0),
    structural_parameters = (kh=kh, ktheta=kθ, m=m, Stheta=Sθ, Itheta=Iθ),
    additional_parameters = (U=U, rho=ρ, c=c)
)

using DifferentialEquations

f = ODEFunction(model)

# initial guess for state variables
x0 = zeros(number_of_states(model))

# steady state problem
prob = SteadyStateProblem(f, x0, p)

# steady state solution
x_ss = solve(prob, SSRootfind())

K, M = linearize(model, x_ss, p)

λ, U, V = get_eigen(model, K, M)

using DifferentialEquations

# construct ODE function
f = ODEFunction(model)

# non-zero plunge degree of freedom
x0 = assemble_states(model;
    aerodynamic_states = (;lambda=zeros(4)),
    structural_states = (;h=0.5, theta=0, hdot=0, thetadot=0))

# simulate for 100 seconds
tspan = (0.0, 100.0)

# assemble problem
prob = DifferentialEquations.ODEProblem(f, x0, tspan, p)

# solve ODE problem
sol = DifferentialEquations.solve(prob)

using Plots
pyplot()

plot(sol,
    vars = [5,6,7,8],
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

# animation parameters
a = -1/5
b = 0.5

# create animation
@gif for t in range(tspan[1], tspan[2], length=200)
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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

