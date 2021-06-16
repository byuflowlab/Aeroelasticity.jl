using StaticArrays, LinearAlgebra

# dimensionless parameters
a = -1/5 # center of mass normalized location
e = -1/10 # reference point normalized location
μ = 20 # = m/(ρ*pi*b^2)
r2 = 6/25 # = Ip/(m*b^2)
σ = 2/5 # = ωh/ωθ
V = range(1, 3, length=100) # = U/(b*ωθ)

# chosen dimensional parameters
b = 2
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

# eigenvalue storage
λ = zeros(ComplexF64, 4, length(V))

for i = 1:length(V)

    M = @SMatrix [1 0 0 0;       #dh/b
              0 1 0 0;       #dθ
              0 0 m m*b*xθ;  #dhdot/b
              0 0 m*b*xθ Ip] #dθdot

    # # mass matrix
    # M = @SMatrix [1 0 0 0;
    #      0 1 0 0;
    #      0 0 b^2/U[i]^2 b^2/U[i]^2*xθ;
    #      0 0 b^2/U[i]^2*xθ b^2/U[i]^2*r2]

    # # mass matrix
    # M = @SMatrix [1 0 0 0;
    #      0 1 0 0;
    #      0 0 1 xθ;
    #      0 0 xθ r2]

    # jacobian
    L_θ = 2*pi*ρ*b*U[i]^2
    J = @SMatrix [0 0 1 0; 0 0 0 1; -kh -L_θ 0 0; 0 b*(1/2+a)*L_θ-kθ 0 0]

    # J = @SMatrix [0 0 1 0;
    #      0 0 0 1;
    #      -m*b^2*ωh^2 -2*pi*ρ*b^2*U[i]^2 0 0;
    #      0 2*(a+1/2)*pi*ρ*b^2*U[i]^2 - Ip*ωθ^2 0 0]

     # # jacobian
     # J = @SMatrix [0 0 1 0;
     #      0 0 0 1;
     #      -σ^2/V[i]^2 -2/μ 0 0;
     #      0 2/μ*(a+1/2)-r2/V[i]^2 0 0]

    # solve generalized eigenvalue problem
    λ[:,i] = eigvals(J, M)# .* U[i]/b
end

using Plots
pyplot()

p1 = scatter(V, imag.(λ')/(ωθ),
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
