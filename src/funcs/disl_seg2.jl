# This file contains the line dislocation from Paul Segall's book

export disp_antiplane, stress_antiplane

function disp_antiplane(x1::T, x2::T, y1::T, d1::T, d2::T, s::T=one(T)) where T
    # must abs(d1) > abs(d2)
    -s / 2π * (
        atan((x1 - y1) / (x2 + d1)) - atan((x1 - y1) / (x2 - d1)) -
        atan((x1 - y1) / (x2 + d2)) + atan((x1 - y1) / (x2 - d2))
    )
end

function stress_antiplane(x1::T, x2::T, y1::T, d1::T, d2::T, G::T, s::T=one(T)) where T
    s13 = -s / 2π * (
        (d1 + x2)/((d1 + x2)^2 + (x1 - y1)^2) - (-(d1 - x2))/((d1 - x2)^2 + (x1 - y1)^2) -
        (d2 + x2)/((d2 + x2)^2 + (x1 - y1)^2) + (-(d2 - x2))/((d2 - x2)^2 + (x1 - y1)^2)
    )
    s23 = -s / 2π * (
        (-(x1 - y1))/((d1 + x2)^2 + (x1 - y1)^2) - (-(x1 - y1))/((d1 - x2)^2 + (x1 - y1)^2) -
        (-(x1 - y1))/((d2 + x2)^2 + (x1 - y1)^2) + (-(x1 - y1))/((d2 - x2)^2 + (x1 - y1)^2)
    )
    return s13, s23
end

function half_infinite_inplane_disp(x1::T, x2::T, ξ1::T, ξ2::T, s1::T, s2::T, ν::T) where T
    # θ1 = atan((x1 - ξ1) / (x2 - ξ2))
    # θ2 = atan((x1 - ξ1) / (x2 + ξ2))
    # r₁ = hypot(x1 - ξ1, x2 - ξ2)
    # r₂ = hypot(x1 - ξ1, x2 + ξ2)
    # r₁² = (x1 - ξ1)^2 + (x2 - ξ2)^2
    # r₂² = (x1 - ξ1)^2 + (x2 + ξ2)^2
    # r₁⁴ = r₁² ^ 2
    # r₂⁴ = r₂² ^ 2
    #
    # u1 = (
    #     s1/π/(1-ν) * ( (1-ν)/2 * (θ2-θ1) + (x1-ξ1)*(x2-ξ2)/4r₁² - (x1-ξ1)*(x2+(3-4ν)*ξ2)/4r₂² + ξ2*x2*x1*(x2+ξ2)/r₂⁴ ) +
    #     s2/π/(1-ν) * ( (1-2ν)/4 * log(r₂/r₁) - (x2-ξ2)^2/4r₁² + (x2^2+ξ2^2-4(1-ν)*ξ2(x2+ξ2))/4r₂² )
    # )
end

function half_infinite_inplane_stress(x1::T, x2::T, ξ1::T, ξ2::T, s1::T, s2::T, μ::T, ν::T) where T
    # r₁² = (x1 - ξ1)^2 + (x2 - ξ2)^2
    # r₂² = (x1 - ξ1)^2 + (x2 + ξ2)^2
end
