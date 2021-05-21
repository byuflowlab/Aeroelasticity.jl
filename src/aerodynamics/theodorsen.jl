using SpecialFunctions

struct Theodorsen <: AerodynamicModel end

function theodorsen(k)
    h0 = hankelh2(0, k)
    h1 = hankelh2(1, k)
    return h1/(h1 + 1im*h0)
end

angle_of_attack(model::Theodorsen, k) = Ck*(θ + hdot/U + b/U*(1/2-a)*θdot)

lift(model::Theodorsen, a, b, ρ, U, θ, θdot, θddot, h, hdot, hddot, Ck) = 2*pi*ρ*U*b*Ck*(hdot + U*θ + b*(1/2-a)*θdot) + pi*ρ*b^2*(hddot + U*θdot - b*a*θddot)

moment(model::Theodorsen, a, b, ρ, U, θdot, θddot, hddot) = -pi*ρ*b^3*(1/2*hddot + U*θdot + b*(1/8 - a/2)*θddot)
