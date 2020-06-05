# This file contains the line dislocation from Paul Segall's book

export disp_antiplane_seg2, stress_antiplane_seg2
export disp_inplane_seg2, stress_inplane_seg2

function disp_antiplane_seg2(x1::T, x2::T, y1::T, d1::T, d2::T, s::T=one(T)) where T
    @assert abs(d1) > abs(d2) "d₁ denotes a deeper dislocation than d₂"
    -s / 2π * (
        atan((x1 - y1) / (x2 + d1)) - atan((x1 - y1) / (x2 - d1)) -
        atan((x1 - y1) / (x2 + d2)) + atan((x1 - y1) / (x2 - d2))
    )
end

function stress_antiplane_seg2(x1::T, x2::T, y1::T, d1::T, d2::T, G::T, s::T=one(T)) where T
    s13 = -s / 2π * G * (
        (d1 + x2)/((d1 + x2)^2 + (x1 - y1)^2) - (-(d1 - x2))/((d1 - x2)^2 + (x1 - y1)^2) -
        (d2 + x2)/((d2 + x2)^2 + (x1 - y1)^2) + (-(d2 - x2))/((d2 - x2)^2 + (x1 - y1)^2)
    )
    s23 = -s / 2π * G * (
        (-(x1 - y1))/((d1 + x2)^2 + (x1 - y1)^2) - (-(x1 - y1))/((d1 - x2)^2 + (x1 - y1)^2) -
        (-(x1 - y1))/((d2 + x2)^2 + (x1 - y1)^2) + (-(x1 - y1))/((d2 - x2)^2 + (x1 - y1)^2)
    )
    return s13, s23
end

function _disp_half_infinite_line_inplane(x1::T, x2::T, ξ1::T, ξ2::T, s1::T, s2::T, ν::T) where T
    θ1 = atan((x1 - ξ1) / (x2 - ξ2))
    θ2 = atan((x1 - ξ1) / (x2 + ξ2))
    r₁ = hypot(x1 - ξ1, x2 - ξ2)
    r₂ = hypot(x1 - ξ1, x2 + ξ2)
    r₁² = r₁ ^ 2
    r₂² = r₂ ^ 2
    r₁⁴ = r₁² ^ 2
    r₂⁴ = r₂² ^ 2

    u1 = (
        -s1/π/(1-ν) * ( (1-ν)/2*(θ2-θ1) + (x1-ξ1)*(x2-ξ2)/4r₁² - (x1-ξ1)*(x2+(3-4ν)*ξ2)/4r₂² + ξ2*x2*(x1-ξ1)*(x2+ξ2)/r₂⁴ ) +
        s2/π/(1-ν) * ( (1-2ν)/4*log(r₂/r₁) - (x2-ξ2)^2/4r₁² + (x2^2+ξ2^2-4(1-ν)*ξ2*(x2+ξ2))/4r₂² + x2*ξ2*(x2+ξ2)^2/r₂⁴ )
    )
    u2 = (
        -s1/π/(1-ν) * ( (1-2ν)/4*log(r₂/r₁) + (x2-ξ2)^2/4r₁² - ((x2+ξ2)^2-2*ξ2^2-2(1-2ν)*ξ2*(x2+ξ2))/4r₂² + x2*ξ2*(x2+ξ2)^2/r₂⁴ ) +
        s2/π/(1-ν) * ( (1-ν)/2*(θ1-θ2) + (x1-ξ1)*(x2-ξ2)/4r₁² - (x1-ξ1)*(x2+(3-4ν)*ξ2)/4r₂² - ξ2*x2*(x1-ξ1)*(x2+ξ2)/r₂⁴ )
    )
    return u1, u2
end

function disp_inplane_seg2(x1::T, x2::T, ξ1::T, ξ2::T, ζ1::T, ζ2::T, s1::T, s2::T, ν::T) where T
    @assert abs(ξ2) > abs(ζ2) "ξ₂ denotes a deeper dislocation than ζ₂"
    u1ξ, u2ξ = _disp_half_infinite_line_inplane(x1, x2, ξ1, ξ2, -s1, -s2, ν)
    u1ζ, u2ζ = _disp_half_infinite_line_inplane(x1, x2, ζ1, ζ2, s1, s2, ν)
    return u1ξ+u1ζ, u2ξ+u2ζ
end

function _stress_half_infinite_line_inplane(x1::T, x2::T, ξ1::T, ξ2::T, s1::T, s2::T, μ::T, ν::T) where T
    r₁ = hypot(x1 - ξ1, x2 - ξ2)
    r₂ = hypot(x1 - ξ1, x2 + ξ2)
    r₁² = r₁ ^ 2
    r₂² = r₂ ^ 2
    r₁⁴ = r₁² ^ 2
    r₂⁴ = r₂² ^ 2
    r₁⁶ = r₁² * r₁⁴
    r₂⁶ = r₂² * r₂⁴

    σ11 = (
        μ*s2/2π/(1-ν) * ( (x1-ξ1)*((x2-ξ2)^2-(x1-ξ1)^2)/r₁⁴ - (x1-ξ1)*((x2+ξ2)^2-(x1-ξ1)^2)/r₂⁴ + 4ξ2*(x1-ξ1)/r₂⁶*((2ξ2-x2)*(x2+ξ2)^2+(3x2+2ξ2)*(x1-ξ1)^2) ) +
        μ*s1/2π/(1-ν) * ( (x2-ξ2)*((x2-ξ2)^2+3(x1-ξ1)^2)/r₁⁴ - (x2+ξ2)*((x2+ξ2)^2+3(x1-ξ1)^2)/r₂⁴ + 2ξ2/r₂⁶*(6x2*(x2+ξ2)*(x1-ξ1)^2-(x2-ξ2)*(x2+ξ2)^3-(x1-ξ1)^4) )
    )

    σ22 = (
        -μ*s2/2π/(1-ν) * ( (x1-ξ1)*(3(x2-ξ2)^2+(x1-ξ1)^2)/r₁⁴ - (x1-ξ1)*(3(x2+ξ2)^2+(x1-ξ1)^2)/r₂⁴ - 4ξ2*x2*(x1-ξ1)/r₂⁶*(3(x2+ξ2)^2-(x1-ξ1)^2) ) +
        μ*s1/2π/(1-ν) * ( (x2-ξ2)*((x2-ξ2)^2-(x1-ξ1)^2)/r₁⁴ - (x2+ξ2)*((x2+ξ2)^2-(x1-ξ1)^2)/r₂⁴ - 2ξ2/r₂⁶*(6x2*(x2+ξ2)*(x1-ξ1)^2-(3x2+ξ2)*(x2+ξ2)^3+(x1-ξ1)^4)  )
    )

    σ12 = (
        μ*s2/2π/(1-ν) * ( (x2-ξ2)*((x2-ξ2)^2-(x1-ξ1)^2)/r₁⁴ - (x2+ξ2)*((x2+ξ2)^2-(x1-ξ1)^2)/r₂⁴ + 2ξ2/r₂⁶*(6x2*(x2+ξ2)*(x1-ξ1)^2-(x1-ξ1)^4+(ξ2-x2)*(x2+ξ2)^3) ) +
        μ*s1/2π/(1-ν) * ( (x1-ξ1)*((x2-ξ2)^2-(x1-ξ1)^2)/r₁⁴ - (x1-ξ1)*((x2+ξ2)^2-(x1-ξ1)^2)/r₂⁴ + 4ξ2*x2*(x1-ξ1)/r₂⁶*(3(x2+ξ2)^2-(x1-ξ1)^2) )
    )

    return σ11, σ22, σ12
end

function stress_inplane_seg2(x1::T, x2::T, ξ1::T, ξ2::T, ζ1::T, ζ2::T, s1::T, s2::T, G::T, ν::T) where T
    @assert abs(ξ2) > abs(ζ2) "ξ₂ denotes a deeper dislocation than ζ₂"
    σ11ξ, σ22ξ, σ12ξ = _stress_half_infinite_line_inplane(x1, x2, ξ1, ξ2, -s1, -s2, G, ν)
    σ11ζ, σ22ζ, σ12ζ = _stress_half_infinite_line_inplane(x1, x2, ζ1, ζ2, s1, s2, G, ν)
    return σ11ξ+σ11ζ, σ22ξ+σ22ζ, σ12ξ+σ12ζ
end
