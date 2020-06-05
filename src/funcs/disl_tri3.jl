# Reference journal article:
# Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
# artefact-free solution.
# Submitted to Geophysical Journal International
#
# Copyright (c) 2014 Mehdi Nikkhoo
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the
# following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# I appreciate any comments or bug reports.
#
# Mehdi Nikkhoo
# created: 2013.1.24
# Last modified: 2014.7.30
#
# VolcanoTectonics Research Group
# Section 2.1, Physics of Earthquakes and Volcanoes
# Department 2, Physics of the Earth
# Helmholtz Centre Potsdam
# German Research Centre for Geosciences (GFZ)
#
# email:
# mehdi.nikkhoo@gfz-potsdam.de
# mehdi.nikkhoo@gmail.com
# -------------------------------------------------------
# Translated by Pengcheng Shi (shipengcheng1230@gmail.com) 08/2018

export disp_tri3_fs, disp_tri3_hs, strain_tri3_fs, strain_tri3_hs, stress_tri3_fs, stress_tri3_hs

const _ey = [0, 1, 0]
const _ez = [0, 0, 1]

"""
    disp_tri3_hs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T) where {T, V}

Compute displacement risen from triangular dislocation in elastic *halfspace*.
    Please see [original version (in supporting information)](https://academic.oup.com/gji/article/201/2/1119/572006#86405752)
    for details, especially the **coordinate system** used here.

## Arguments
- `X`, `Y`, `Z`: observational coordinates
- `P1`, `P2`, `P3`: three triangular vertices coordinates respectively
- `Ss`, `Ds`, `Ts`: triangular dislocation vector, Strike-slip, Dip-slip, Tensile-slip respectively
- `nu`: poisson ratio

## Output
By order: ``u_x``, ``u_y``, ``u_z``
"""
function disp_tri3_hs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T) where {T, V}
    @assert (Z ≤ zero(T) && P1[3] ≤ zero(T) && P2[3] ≤ zero(T) && P3[3] ≤ zero(T))  "Half-space solution: Z coordinates must be negative!"
    A = transform_matrix(P1, P2, P3)
    ueMS, unMS, uvMS = _disp_tri3_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A')
    ueFSC, unFSC, uvFSC = _td_disp_harmonic_func(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A)
    P1[3] *= -one(T)
    P2[3] *= -one(T)
    P3[3] *= -one(T)
    allatsurface = P1[3] == zero(T) && P2[3] == zero(T) && P3[3] == zero(T)
    if !allatsurface
        A *= -one(T)
        A[3,1] *= -one(T)
        A[3,2] *= -one(T)
        A[3,3] *= -one(T)
    end
    ueIS, unIS, uvIS = _disp_tri3_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A')
    uvIS = ifelse(allatsurface, -uvIS, uvIS)
    ue = ueMS + ueIS + ueFSC
    un = unMS + unIS + unFSC
    uv = uvMS + uvIS + uvFSC
    allatsurface && begin ue *= -one(T); un *= -one(T); uv *= -one(T) end
    P1[3] *= -one(T)
    P2[3] *= -one(T)
    P3[3] *= -one(T)
    return ue, un, uv
end

"""
    disp_tri3_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T) where {T, V}

Compute displacement risen from triangular dislocation in elastic *fullspace*.
    Please see [original version (in supporting information)](https://academic.oup.com/gji/article/201/2/1119/572006#86405752)
    for details, especially the **coordinate system** used here.

## Arguments
The same as [`disp_tri3_hs`](@ref)
"""
function disp_tri3_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T) where {T, V}
    A = transform_matrix(P1, P2, P3)
    _disp_tri3_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu, A')
end

@inline function transform_matrix(P1::V, P2::V, P3::V) where V<:AbstractVector{T} where T
    Vnorm = cross(P2 - P1, P3 - P1)
    normalize!(Vnorm)
    Vstrike = cross(_ez, Vnorm)

    if norm(Vstrike) == zero(T)
        Vstrike = _ey * Vnorm[3]
        if P1[3] > zero(T)
            Vstrike *= -one(T)
        end
    end
    normalize!(Vstrike)
    Vdip = Vnorm × Vstrike
    hcat(Vnorm, Vstrike, Vdip)
end

@inline function _disp_tri3_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T, At::M) where {T, V, M}
    bx, by, bz = Ts, Ss, Ds
    p1, p2, p3 = [zeros(T, 3) for _ in 1: 3]

    x, y, z = coord_trans(X - P2[1], Y - P2[2], Z - P2[3], At)
    p1[1], p1[2], p1[3] = coord_trans(P1[1] - P2[1], P1[2] - P2[2], P1[3] - P2[3], At)
    p3[1], p3[2], p3[3] = coord_trans(P3[1] - P2[1], P3[2] - P2[2], P3[3] - P2[3], At)

    e12 = normalize(p2 - p1)
    e13 = normalize(p3 - p1)
    e23 = normalize(p3 - p2)

    A = acos(dot(e12, e13))
    B = acos(dot(-e12, e23))
    C = acos(dot(e23, e13))

    Trimode = trimodefinder(y, z, x, p1, p2, p3)

    if Trimode == 1
        u1Tp, v1Tp, w1Tp = TDSetupD(x, y, z, A, bx, by, bz, nu, p1, -e13)
        u2Tp, v2Tp, w2Tp = TDSetupD(x, y, z, B, bx, by, bz, nu, p2, e12)
        u3Tp, v3Tp, w3Tp = TDSetupD(x, y, z, C, bx, by, bz, nu, p3, e23)
        u = u1Tp + u2Tp + u3Tp
        v = v1Tp + v2Tp + v3Tp
        w = w1Tp + w2Tp + w3Tp
    end

    if Trimode == -1
        u1Tn, v1Tn, w1Tn = TDSetupD(x, y, z, A, bx, by, bz, nu, p1, e13)
        u2Tn, v2Tn, w2Tn = TDSetupD(x, y, z, B, bx, by, bz, nu, p2, -e12)
        u3Tn, v3Tn, w3Tn = TDSetupD(x, y, z, C, bx, by, bz, nu, p3, -e23)
        u = u1Tn + u2Tn + u3Tn
        v = v1Tn + v2Tn + v3Tn
        w = w1Tn + w2Tn + w3Tn
    end

    if Trimode == 0
        u = NaN
        v = NaN
        w = NaN
    end

    a1, a2, a3 = -x, p1[2] - y, p1[3] - z
    b1, b2, b3 = -x, -y, -z
    c1, c2, c3 = -x, p3[2] - y, p3[3] - z
    na = hypot(a1, a2, a3)
    nb = hypot(b1, b2, b3)
    nc = hypot(c1, c2, c3)

    Fiy = a1 * (b2 * c3 - b3 * c2)- a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)
    Fix = na * nb * nc + (a1 * b1 + a2 * b2 + a3 * b3) * nc + (a1 * c1 + a2 * c2 + a3 * c3) * nb + (b1 * c1 + b2 * c2 + b3 * c3) * na
    Fi = -2 * atan(Fiy, Fix) / 4 / π
    u += bx * Fi
    v += by * Fi
    w += bz * Fi
    ue, un, uv = coord_trans(u, v, w, At')
end

@inline function _td_disp_harmonic_func(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, nu::T, A::M) where {T, V, M}
    bx, by, bz = Ts, Ss, Ds
    bX, bY, bZ = coord_trans(bx, by, bz, A)
    u1, v1, w1 = AngSetupFSC(X, Y, Z, bX, bY, bZ, P1, P2, nu)
    u2, v2, w2 = AngSetupFSC(X, Y, Z, bX, bY, bZ, P2, P3, nu)
    u3, v3, w3 = AngSetupFSC(X, Y, Z, bX, bY, bZ, P3, P1, nu)
    return u1 + u2 + u3, v1 + v2 + v3, w1 + w2 + w3
end

@inline function AngSetupFSC(X::T, Y::T, Z::T, bX::T, bY::T, bZ::T, PA::V, PB::V, nu::T) where {T, V}
    sv1, sv2, sv3 = PB[1] - PA[1], PB[2] - PA[2], PB[3] - PA[3]
    svr = hypot(sv1, sv2, sv3)
    beta = acos(-sv3 / svr)
    if beta ≈ zero(T) || beta ≈ π
        return zero(T), zero(T), zero(T)
    else
        ey1 = [sv1, sv2, zero(T)]
        normalize!(ey1)
        ey3 = -_ez
        ey2 = cross(ey3, ey1)
        A = hcat(ey1, ey2, ey3)
        y1A, y2A, y3A = coord_trans(X - PA[1], Y - PA[2], Z - PA[3], A)
        y1AB, y2AB, y3AB = coord_trans(sv1, sv2, sv3, A)
        y1B = y1A - y1AB
        y2B = y2A - y2AB
        y3B = y3A - y3AB
        b1, b2, b3 = coord_trans(bX, bY, bZ, A)
        angle = ifelse(beta * y1A ≥ zero(T), -π + beta, beta)
        v1A, v2A, v3A = AngDisDispFSC(y1A, y2A, y3A, angle, b1, b2, b3, nu, -PA[3])
        v1B, v2B, v3B = AngDisDispFSC(y1B, y2B, y3B, angle, b1, b2, b3, nu, -PB[3])
        v1 = v1B - v1A
        v2 = v2B - v2A
        v3 = v3B - v3A
        ue, un, uv = coord_trans(v1, v2, v3, A')
        return ue, un, uv
    end
end

@inline function AngDisDispFSC(y1::T, y2::T, y3::T, beta::T, b1::T, b2::T, b3::T, nu::T, a::T) where T
    sinB, cosB = sincos(beta)
    cotB = cot(beta)
    y3b = y3 + 2a
    z1b = y1 * cosB + y3b * sinB
    z3b = -y1 * sinB + y3b * cosB
    rb = hypot(y1, y2, y3b)
    Fib = 2 * atan(-y2 / (-(rb + y3b) * cot(beta/2) + y1))

    v1cb1 = (b1 / 4 / pi / (1 - nu) * (-2 * (1 - nu) * (1 - 2 * nu) * Fib * cotB ^ 2 + (1 - 2 * nu) * y2 /
        (rb + y3b) * ((1 - 2 * nu - a / rb) * cotB - y1 / (rb + y3b) * (nu + a / rb)) + (1 - 2 * nu) *
        y2 * cosB * cotB / (rb + z3b) * (cosB + a / rb) + a * y2 * (y3b - a) * cotB / rb ^ 3 + y2 *
        (y3b - a) / (rb * (rb + y3b)) * (-(1 - 2 * nu) * cotB + y1 / (rb + y3b) * (2 * nu + a / rb) +
        a * y1 / rb ^ 2) +y2 * (y3b - a) / (rb * (rb + z3b)) * (cosB / (rb + z3b) * ((rb *
        cosB + y3b) * ((1 - 2 * nu) * cosB - a / rb) * cotB + 2 * (1 - nu) * (rb * sinB - y1) * cosB) -
        a * y3b * cosB * cotB / rb ^ 2)))

    v2cb1 = ((b1 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * ((2 * (1 - nu) * cotB ^ 2 - nu) * log(rb + y3b) - (2 *
        (1 - nu) * cotB ^ 2 + 1 - 2 * nu)* cosB * log(rb + z3b)) -(1 - 2 * nu) / (rb + y3b) * (y1 *
        cotB * (1 - 2 * nu - a / rb) + nu * y3b - a + y2 ^ 2 / (rb + y3b) * (nu + a / rb)) -(1 - 2 *
        nu) * z1b * cotB / (rb + z3b) * (cosB + a / rb) - a * y1 * (y3b - a) * cotB / rb ^ 3 +
        (y3b - a) / (rb + y3b) * (-2 * nu + 1 / rb * ((1 - 2 * nu) * y1 * cotB - a) + y2 ^ 2 / (rb *
        (rb + y3b)) * (2 * nu + a / rb) + a * y2 ^ 2 / rb ^ 3) +(y3b - a) / (rb + z3b) * (cosB ^ 2 -
        1 / rb * ((1 - 2 * nu) * z1b * cotB + a * cosB) + a * y3b * z1b * cotB / rb ^ 3 - 1 / (rb *
        (rb + z3b)) * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB / rb * (rb * cosB + y3b))))))

    v3cb1 = (b1 / 4 / pi / (1 - nu) * (2 * (1 - nu) * (((1 - 2 * nu) * Fib * cotB) + (y2 / (rb + y3b) * (2 *
        nu + a / rb)) -(y2 * cosB / (rb + z3b) * (cosB + a / rb))) +y2 * (y3b - a) / rb * (2 *
        nu / (rb + y3b) + a / rb ^ 2) +y2 * (y3b - a) * cosB / (rb * (rb + z3b)) * (1 - 2 * nu -
        (rb * cosB + y3b) / (rb + z3b) * (cosB + a / rb) - a * y3b / rb ^ 2)))

    v1cb2 = (b2 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * ((2 * (1 - nu) * cotB ^ 2 + nu) * log(rb + y3b) - (2 *
        (1 - nu) * cotB ^ 2 + 1)* cosB * log(rb + z3b)) +(1 - 2 * nu) / (rb + y3b) * (-(1 - 2 * nu) *
        y1 * cotB + nu * y3b - a + a * y1 * cotB / rb + y1 ^ 2 / (rb + y3b) * (nu + a / rb)) -(1 - 2 *
        nu)* cotB / (rb + z3b) * (z1b * cosB - a * (rb * sinB - y1) / (rb * cosB)) - a * y1 *
        (y3b - a) * cotB / rb ^ 3 + (y3b - a) / (rb + y3b) * (2 * nu + 1 / rb * ((1 - 2 * nu) * y1 *
        cotB + a) -y1 ^ 2 / (rb * (rb + y3b)) * (2 * nu + a / rb) - a * y1 ^ 2 / rb ^ 3) +(y3b - a) *
        cotB / (rb + z3b) * (-cosB * sinB + a * y1 * y3b / (rb ^ 3 * cosB) + (rb * sinB - y1) /
        rb * (2 * (1 - nu) * cosB - (rb * cosB + y3b) / (rb + z3b) * (1 + a / (rb * cosB))))))

    v2cb2 = (b2 / 4 / pi / (1 - nu) * (2 * (1 - nu) * (1 - 2 * nu) * Fib * cotB ^ 2 + (1 - 2 * nu) * y2 /
        (rb + y3b) * (-(1 - 2 * nu - a / rb) * cotB + y1 / (rb + y3b) * (nu + a / rb)) - (1 - 2 * nu) *
        y2 * cotB / (rb + z3b) * (1 + a / (rb * cosB)) - a * y2 * (y3b - a) * cotB / rb ^ 3 + y2 *
        (y3b - a) / (rb * (rb + y3b)) * ((1 - 2 * nu) * cotB - 2 * nu * y1 / (rb + y3b) - a * y1 / rb *
        (1 / rb + 1 / (rb + y3b))) +y2 * (y3b - a) * cotB / (rb * (rb + z3b)) * (-2 * (1 - nu) *
        cosB + (rb * cosB + y3b) / (rb + z3b) * (1 + a / (rb * cosB)) + a * y3b / (rb ^ 2 * cosB))))

    v3cb2 = (b2 / 4 / pi / (1 - nu) * (-2 * (1 - nu) * (1 - 2 * nu) * cotB * (log(rb + y3b) - cosB *
        log(rb + z3b)) -2 * (1 - nu) * y1 / (rb + y3b) * (2 * nu + a / rb) + 2 * (1 - nu) * z1b / (rb +
        z3b) * (cosB + a / rb) + (y3b - a) / rb * ((1 - 2 * nu) * cotB - 2 * nu * y1 / (rb + y3b) - a *
        y1 / rb ^ 2) -(y3b - a) / (rb + z3b) * (cosB * sinB + (rb * cosB + y3b) * cotB / rb *
        (2 * (1 - nu) * cosB - (rb * cosB + y3b) / (rb + z3b)) + a / rb * (sinB - y3b * z1b /
        rb ^ 2 - z1b * (rb * cosB + y3b) / (rb * (rb + z3b))))))

    v1cb3 = (b3 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * (y2 / (rb + y3b) * (1 + a / rb) - y2 * cosB / (rb +
        z3b) * (cosB + a / rb)) -y2 * (y3b - a) / rb * (a / rb ^ 2 + 1 / (rb + y3b)) + y2 *
        (y3b - a) * cosB / (rb * (rb + z3b)) * ((rb * cosB + y3b) / (rb + z3b) * (cosB + a /
        rb) +a * y3b / rb ^ 2)))

    v2cb3 = (b3 / 4 / pi / (1 - nu) * ((1 - 2 * nu) * (-sinB * log(rb + z3b) - y1 / (rb + y3b) * (1 + a /
        rb) +z1b / (rb + z3b) * (cosB + a / rb)) +y1 * (y3b - a) / rb * (a / rb ^ 2 + 1 / (rb +
        y3b)) -(y3b - a) / (rb + z3b) * (sinB * (cosB - a / rb) + z1b / rb * (1 + a * y3b /
        rb ^ 2) -1 / (rb * (rb + z3b)) * (y2 ^ 2 * cosB * sinB - a * z1b / rb * (rb * cosB + y3b)))))

    v3cb3 = (b3 / 4 / pi / (1 - nu) * (2 * (1 - nu) * Fib + 2 * (1 - nu) * (y2 * sinB / (rb + z3b) * (cosB +
        a / rb)) +y2 * (y3b - a) * sinB / (rb * (rb + z3b)) * (1 + (rb * cosB + y3b) / (rb +
        z3b) * (cosB + a / rb) + a * y3b / rb ^ 2)))

    return v1cb1 + v1cb2 + v1cb3, v2cb1 + v2cb2 + v2cb3, v3cb1 + v3cb2 + v3cb3
end

@inline function trimodefinder(x::T, y::T, z::T, p1::V, p2::V, p3::V) where {T, V}
    a = ((p2[3] - p3[3]) * (x - p3[2]) + (p3[2] - p2[2]) * (y - p3[3])) /
        ((p2[3] - p3[3]) * (p1[2] - p3[2]) + (p3[2] - p2[2]) * (p1[3] - p3[3]))
    b = ((p3[3] - p1[3]) * (x - p3[2]) + (p1[2] - p3[2]) * (y - p3[3])) /
        ((p2[3] - p3[3]) * (p1[2] - p3[2]) + (p3[2] - p2[2]) * (p1[3] - p3[3]))
    c = 1 - a - b
    trimode = 1
    if ((a ≤ 0) && (b > c) && (c > a)) || ((b ≤ 0) && (c > a) && (a > b)) || ((c ≤ 0) && (a > b) && (b > c))
        trimode = -1
    end
    if ((a == 0) && (b ≥ 0) && (c ≥ 0)) || ((a ≥ 0) && (b == 0) && (c ≥ 0)) || ((a ≥ 0) && (b ≥ 0) && (c == 0))
        trimode = 0
    end
    if (trimode == 0) && (z ≠ 0)
        trimode = 1
    end
    return trimode
end

@inline function TDSetupD(x::T, y::T, z::T, alpha::T, bx::T, by::T, bz::T, nu::T, TriVertex::V, SideVec::V) where {T, V}
    y1 = SideVec[3] * (y - TriVertex[2]) - SideVec[2] * (z - TriVertex[3])
    z1 = SideVec[2] * (y - TriVertex[2]) + SideVec[3] * (z - TriVertex[3])
    by1 = SideVec[3] * by - SideVec[2] * bz
    bz1 = SideVec[2] * by + SideVec[3] * bz
    u, v0, w0 = AngDisDisp(x, y1, z1, -π + alpha, bx, by1, bz1, nu)
    v = SideVec[3] * v0 + SideVec[2] * w0
    w = -SideVec[2] * v0 + SideVec[3] * w0
    return u, v, w
end

@inline function _TDstrain_HarFunc(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T, At::M) where {T, V, M}
    bx, by, bz = Ts, Ss, Ds
    bX, bY, bZ = coord_trans(bx, by, bz, At)
    exx1, eyy1, ezz1, exy1, exz1, eyz1 = AngSetupFSC_S(X, Y, Z, bX, bY, bZ, P1, P2, λ, μ)
    exx2, eyy2, ezz2, exy2, exz2, eyz2 = AngSetupFSC_S(X, Y, Z, bX, bY, bZ, P2, P3, λ, μ)
    exx3, eyy3, ezz3, exy3, exz3, eyz3 = AngSetupFSC_S(X, Y, Z, bX, bY, bZ, P3, P1, λ, μ)
    return exx1 + exx2 + exx3, eyy1 + eyy2 + eyy3, ezz1 + ezz2 + ezz3, exy1 + exy2 + exy3, exz1 + exz2 + exz3, eyz1 + eyz2 + eyz3
end

@inline function AngSetupFSC_S(X::T, Y::T, Z::T, bX::T, bY::T, bZ::T, PA::V, PB::V, λ::T, μ::T) where {T, V}
    nu = λ / (μ + λ) / 2
    sv1, sv2, sv3 = PB[1] - PA[1], PB[2] - PA[2], PB[3] - PA[3]
    svr = hypot(sv1, sv2, sv3)
    beta = acos(-sv3 / svr)
    if beta ≈ zero(T) || beta ≈ π
        return zero(T), zero(T), zero(T), zero(T), zero(T), zero(T)
    else
        ey1 = [sv1, sv2, zero(T)]
        normalize!(ey1)
        ey3 = -_ez
        ey2 = cross(ey3, ey1)
        A = hcat(ey1, ey2, ey3)
        y1A, y2A, y3A = coord_trans(X - PA[1], Y - PA[2], Z - PA[3], A)
        y1AB, y2AB, y3AB = coord_trans(sv1, sv2, sv3, A)
        y1B = y1A - y1AB
        y2B = y2A - y2AB
        y3B = y3A - y3AB
        b1, b2, b3 = coord_trans(bX, bY, bZ, A)
        if beta * y1A ≥ zero(T)
            v11A, v22A, v33A, v12A, v13A, v23A = AngDisStrainFSC(-y1A, -y2A, y3A, π - beta, -b1, -b2, b3, nu, -PA[3])
            v13A *= -one(T)
            v23A *= -one(T)
            v11B, v22B, v33B, v12B, v13B, v23B = AngDisStrainFSC(-y1B, -y2B, y3B, π - beta, -b1, -b2, b3, nu, -PB[3])
            v13B *= -one(T)
            v23B *= -one(T)
        else
            v11A, v22A, v33A, v12A, v13A, v23A = AngDisStrainFSC(y1A, y2A, y3A, beta, b1, b2, b3, nu, -PA[3])
            v11B, v22B, v33B, v12B, v13B, v23B = AngDisStrainFSC(y1B, y2B, y3B, beta, b1, b2, b3, nu, -PB[3])
        end
        v11 = v11B - v11A
        v22 = v22B - v22A
        v33 = v33B - v33A
        v12 = v12B - v12A
        v13 = v13B - v13A
        v23 = v23B - v23A
        Exx, Eyy, Ezz, Exy, Exz, Eyz = TensTrans(v11, v22, v33, v12, v13, v23, A')
        return Exx, Eyy, Ezz, Exy, Exz, Eyz
    end
end

@inline function AngDisStrainFSC(y1::T, y2::T, y3::T, beta::T, b1::T, b2::T, b3::T, nu::T, a::T) where T
    sinB, cosB = sincos(beta)
    cotB = cot(beta)
    y3b = y3 + 2a
    z1b = y1 * cosB + y3b * sinB
    z3b = -y1 * sinB + y3b * cosB
    rb2 = y1^2 + y2^2 + y3b^2
    rb = sqrt(rb2)
    W1 = rb * cosB + y3b
    W2 = cosB + a / rb
    W3 = cosB + y3b / rb
    W4 = nu + a / rb
    W5 = 2nu + a / rb
    W6 = rb + y3b
    W7 = rb + z3b
    W8 = y3 + a
    W9 = 1 + a / rb / cosB
    N1 = 1-2*nu

    rFib_ry2 = z1b / rb / (rb + z3b) - y1 / rb / (rb + y3b) # y2 = x in ADCS
    rFib_ry1 = y2 / rb / (rb + y3b) - cosB * y2 / rb / (rb + z3b) # y1 =y in ADCS
    rFib_ry3 = -sinB * y2 / rb / (rb + z3b) # y3 = z in ADCS

    v11 = (b1 * (1 / 4 * ((-2 + 2 * nu) * N1 * rFib_ry1 * cotB ^ 2 - N1 * y2 / W6 ^ 2 * ((1 - W5) * cotB -
        y1 / W6 * W4) /rb *y1+N1 *y2 /W6 * (a / rb ^ 3 * y1 * cotB - 1 / W6 * W4 + y1 ^ 2 /
        W6 ^ 2 * W4 / rb + y1 ^ 2 / W6 * a / rb ^ 3) -N1 * y2 * cosB * cotB / W7 ^ 2 * W2 * (y1 /
        rb - sinB) -N1 * y2 * cosB * cotB / W7 * a / rb ^ 3 * y1 - 3 * a * y2 * W8 * cotB / rb ^ 5. *
        y1 - y2 * W8 / rb ^ 3 / W6 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) * y1 - y2 * W8 /
        rb2 / W6 ^ 2 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) * y1 + y2 * W8 / rb / W6 *
        (1 / W6 * W5 - y1 ^ 2 / W6 ^ 2 * W5 / rb - y1 ^ 2 / W6 * a / rb ^ 3 + a / rb2 - 2 * a * y1 ^
        2 / rb2 ^ 2) -y2 * W8 / rb ^ 3 / W7 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB +
        (2 - 2 * nu) * (rb * sinB - y1) * cosB) -a * y3b * cosB * cotB / rb2) * y1 - y2 * W8 / rb /
        W7 ^ 2 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + (2 - 2 * nu) * (rb * sinB - y1) *
        cosB) -a * y3b * cosB * cotB / rb2) * (y1 / rb - sinB) + y2 * W8 / rb / W7 * (-cosB /
        W7 ^ 2 * (W1 * (N1 * cosB - a / rb) * cotB + (2 - 2 * nu) * (rb * sinB - y1) * cosB) * (y1 /
        rb - sinB) +cosB / W7 * (1 / rb * cosB * y1 * (N1 * cosB - a / rb) * cotB + W1 * a / rb ^
        3 * y1 * cotB + (2 - 2 * nu) * (1 / rb * sinB * y1 - 1) * cosB) +2 * a * y3b * cosB * cotB /
        rb2 ^ 2 * y1)) /π/(1 - nu)) +
        b2 * (1 / 4 * (N1 * (((2 - 2 * nu) * cotB ^ 2 + nu) / rb * y1 / W6 - ((2 - 2 * nu) * cotB ^ 2 + 1) *
        cosB * (y1 / rb - sinB) / W7) -N1 / W6 ^ 2 * (-N1 * y1 * cotB + nu * y3b - a + a * y1 *
        cotB / rb + y1 ^ 2 / W6 * W4) /rb *y1+N1 /W6 * (-N1 * cotB + a * cotB / rb - a *
        y1 ^ 2 * cotB / rb ^ 3 + 2 * y1 / W6 * W4 - y1 ^ 3 / W6 ^ 2 * W4 / rb - y1 ^ 3 / W6 * a /
        rb ^ 3) +N1 * cotB / W7 ^ 2 * (z1b * cosB - a * (rb * sinB - y1) / rb / cosB) * (y1 /
        rb - sinB) -N1 * cotB / W7 * (cosB ^ 2 - a * (1 / rb * sinB * y1 - 1) / rb / cosB + a *
        (rb * sinB - y1) / rb ^ 3 / cosB * y1) -a * W8 * cotB / rb ^ 3 + 3 * a * y1 ^ 2 * W8 *
        cotB / rb ^ 5 - W8 / W6 ^ 2 * (2 * nu + 1 / rb * (N1 * y1 * cotB + a) - y1 ^ 2 / rb / W6 *
        W5 - a * y1 ^ 2 / rb ^ 3) /rb *y1+W8 /W6 * (-1 / rb ^ 3 * (N1 * y1 * cotB + a) * y1 +
        1 / rb * N1 * cotB - 2 * y1 / rb / W6 * W5 + y1 ^ 3 / rb ^ 3 / W6 * W5 + y1 ^ 3 / rb2 /
        W6 ^ 2 * W5 + y1 ^ 3 / rb2 ^ 2 / W6 * a - 2 * a / rb ^ 3 * y1 + 3 * a * y1 ^ 3 / rb ^ 5) -W8 *
        cotB / W7 ^ 2 * (-cosB * sinB + a * y1 * y3b / rb ^ 3 / cosB + (rb * sinB - y1) / rb *
        ((2 - 2 * nu) * cosB - W1 / W7 * W9)) * (y1 / rb - sinB) + W8 * cotB / W7 * (a * y3b /
        rb ^ 3 / cosB - 3 * a * y1 ^ 2 * y3b / rb ^ 5. / cosB + (1 / rb * sinB * y1 - 1) / rb *
        ((2 - 2 * nu) * cosB - W1 / W7 * W9) - (rb * sinB - y1) / rb ^ 3 * ((2 - 2 * nu) * cosB - W1 /
        W7 * W9) * y1 + (rb * sinB - y1) / rb * (-1 / rb * cosB * y1 / W7 * W9 + W1 / W7 ^ 2 *
        W9 * (y1 / rb - sinB) + W1 / W7 * a / rb ^ 3 / cosB * y1))) /π/(1 - nu)) +
        b3 * (1 / 4 * (N1 * (-y2 / W6 ^ 2 * (1 + a / rb) / rb * y1 - y2 / W6 * a / rb ^ 3 * y1 + y2 *
        cosB / W7 ^ 2 * W2 * (y1 / rb - sinB) + y2 * cosB / W7 * a / rb ^ 3 * y1) +y2 * W8 /
        rb ^ 3 * (a / rb2 + 1 / W6) * y1 - y2 * W8 / rb * (-2 * a / rb2 ^ 2 * y1 - 1 / W6 ^ 2 /
        rb * y1) -y2 * W8 * cosB / rb ^ 3 / W7 * (W1 / W7 * W2 + a * y3b / rb2) * y1 - y2 * W8 *
        cosB / rb / W7 ^ 2 * (W1 / W7 * W2 + a * y3b / rb2) * (y1 / rb - sinB) + y2 * W8 *
        cosB / rb / W7 * (1 / rb * cosB * y1 / W7 * W2 - W1 / W7 ^ 2 * W2 * (y1 / rb - sinB) -
        W1 / W7 * a / rb ^ 3 * y1 - 2 * a * y3b / rb2 ^ 2 * y1)) /π/(1 - nu)))

    v22 = (b1 * (1 / 4 * (N1 * (((2 - 2 * nu) * cotB ^ 2 - nu) / rb * y2 / W6 - ((2 - 2 * nu) * cotB ^ 2 + 1 -
        2 * nu)* cosB / rb * y2 / W7) +N1 / W6 ^ 2 * (y1 * cotB * (1 - W5) + nu * y3b - a + y2 ^
        2 / W6 * W4) /rb *y2-N1 /W6 * (a * y1 * cotB / rb ^ 3 * y2 + 2 * y2 / W6 * W4 - y2 ^
        3 / W6 ^ 2 * W4 / rb - y2 ^ 3 / W6 * a / rb ^ 3) +N1 * z1b * cotB / W7 ^ 2 * W2 / rb *
        y2 + N1 * z1b * cotB / W7 * a / rb ^ 3 * y2 + 3 * a * y2 * W8 * cotB / rb ^ 5. * y1 - W8 /
        W6 ^ 2 * (-2 * nu + 1 / rb * (N1 * y1 * cotB - a) + y2 ^ 2 / rb / W6 * W5 + a * y2 ^ 2 /
        rb ^ 3) /rb *y2+W8 /W6 * (-1 / rb ^ 3 * (N1 * y1 * cotB - a) * y2 + 2 * y2 / rb /
        W6 * W5 - y2 ^ 3 / rb ^ 3 / W6 * W5 - y2 ^ 3 / rb2 / W6 ^ 2 * W5 - y2 ^ 3 / rb2 ^ 2 / W6 *
        a + 2 * a / rb ^ 3 * y2 - 3 * a * y2 ^ 3 / rb ^ 5) -W8 / W7 ^ 2 * (cosB ^ 2 - 1 / rb * (N1 *
        z1b * cotB + a * cosB) +a * y3b * z1b * cotB / rb ^ 3 - 1 / rb / W7 * (y2 ^ 2 * cosB ^ 2 -
        a * z1b * cotB / rb * W1)) /rb *y2+W8 /W7 * (1 / rb ^ 3 * (N1 * z1b * cotB + a *
        cosB) * y2 - 3 * a * y3b * z1b * cotB / rb ^ 5. * y2 + 1 / rb ^ 3 / W7 * (y2 ^ 2 * cosB ^ 2 -
        a * z1b * cotB / rb * W1) * y2 + 1 / rb2 / W7 ^ 2 * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB /
        rb * W1) * y2 - 1 / rb / W7 * (2 * y2 * cosB ^ 2 + a * z1b * cotB / rb ^ 3 * W1 * y2 - a *
        z1b * cotB / rb2 * cosB * y2))) /π/(1 - nu)) +
        b2 * (1 / 4 * ((2 - 2 * nu) * N1 * rFib_ry2 * cotB ^ 2 + N1 / W6 * ((W5 - 1) * cotB + y1 / W6 *
        W4) -N1 * y2 ^ 2 / W6 ^ 2 * ((W5 - 1) * cotB + y1 / W6 * W4) / rb + N1 * y2 / W6 * (-a /
        rb ^ 3 * y2 * cotB - y1 / W6 ^ 2 * W4 / rb * y2 - y2 / W6 * a / rb ^ 3 * y1) -N1 * cotB /
        W7 * W9 + N1 * y2 ^ 2 * cotB / W7 ^ 2 * W9 / rb + N1 * y2 ^ 2 * cotB / W7 * a / rb ^ 3 /
        cosB - a * W8 * cotB / rb ^ 3 + 3 * a * y2 ^ 2 * W8 * cotB / rb ^ 5 + W8 / rb / W6 * (N1 *
        cotB - 2 * nu * y1 / W6 - a * y1 / rb * (1 / rb + 1 / W6)) -y2 ^ 2 * W8 / rb ^ 3 / W6 *
        (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb * (1 / rb + 1 / W6)) - y2 ^ 2 * W8 / rb2 / W6 ^
        2 * (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb * (1 / rb + 1 / W6)) + y2 * W8 / rb / W6 *
        (2 * nu * y1 / W6 ^ 2 / rb * y2 + a * y1 / rb ^ 3 * (1 / rb + 1 / W6) * y2 - a * y1 / rb *
        (-1 / rb ^ 3 * y2 - 1 / W6 ^ 2 / rb * y2)) +W8 * cotB / rb / W7 * ((-2 + 2 * nu) * cosB +
        W1 / W7 * W9 + a * y3b / rb2 / cosB) -y2 ^ 2 * W8 * cotB / rb ^ 3 / W7 * ((-2 + 2 * nu) *
        cosB + W1 / W7 * W9 + a * y3b / rb2 / cosB) -y2 ^ 2 * W8 * cotB / rb2 / W7 ^ 2 * ((-2 +
        2 * nu)* cosB + W1 / W7 * W9 + a * y3b / rb2 / cosB) +y2 * W8 * cotB / rb / W7 * (1 /
        rb * cosB * y2 / W7 * W9 - W1 / W7 ^ 2 * W9 / rb * y2 - W1 / W7 * a / rb ^ 3 / cosB * y2 -
        2 * a * y3b / rb2 ^ 2 / cosB * y2)) /π/(1 - nu)) +
        b3 * (1 / 4 * (N1 * (-sinB / rb * y2 / W7 + y2 / W6 ^ 2 * (1 + a / rb) / rb * y1 + y2 / W6 *
        a / rb ^ 3 * y1 - z1b / W7 ^ 2 * W2 / rb * y2 - z1b / W7 * a / rb ^ 3 * y2) -y2 * W8 /
        rb ^ 3 * (a / rb2 + 1 / W6) * y1 + y1 * W8 / rb * (-2 * a / rb2 ^ 2 * y2 - 1 / W6 ^ 2 /
        rb * y2) +W8 / W7 ^ 2 * (sinB * (cosB - a / rb) + z1b / rb * (1 + a * y3b / rb2) - 1 /
        rb / W7 * (y2 ^ 2 * cosB * sinB - a * z1b / rb * W1)) /rb *y2-W8 /W7 * (sinB * a /
        rb ^ 3 * y2 - z1b / rb ^ 3 * (1 + a * y3b / rb2) * y2 - 2 * z1b / rb ^ 5 * a * y3b * y2 +
        1 / rb ^ 3 / W7 * (y2 ^ 2 * cosB * sinB - a * z1b / rb * W1) * y2 + 1 / rb2 / W7 ^ 2 *
        (y2 ^ 2 * cosB * sinB - a * z1b / rb * W1) * y2 - 1 / rb / W7 * (2 * y2 * cosB * sinB + a *
        z1b / rb ^ 3 * W1 * y2 - a * z1b / rb2 * cosB * y2))) /π/(1 - nu)))

    v33 = (b1 * (1 / 4 * ((2 - 2 * nu) * (N1 * rFib_ry3 * cotB - y2 / W6 ^ 2 * W5 * (y3b / rb + 1) -
        1 / 2 * y2 / W6 * a / rb ^ 3 * 2 * y3b + y2 * cosB / W7 ^ 2 * W2 * W3 + 1 / 2 * y2 * cosB / W7 *
        a / rb ^ 3 * 2 * y3b) +y2 / rb * (2 * nu / W6 + a / rb2) - 1 / 2 * y2 * W8 / rb ^ 3 * (2 *
        nu / W6 + a / rb2)* 2 * y3b + y2 * W8 / rb * (-2 * nu / W6 ^ 2 * (y3b / rb + 1) - a /
        rb2 ^ 2 * 2 * y3b) +y2 * cosB / rb / W7 * (1 - 2 * nu - W1 / W7 * W2 - a * y3b / rb2) -
        1 / 2 * y2 * W8 * cosB / rb ^ 3 / W7 * (1 - 2 * nu - W1 / W7 * W2 - a * y3b / rb2) * 2 *
        y3b - y2 * W8 * cosB / rb / W7 ^ 2 * (1 - 2 * nu - W1 / W7 * W2 - a * y3b / rb2) * W3 + y2 *
        W8 * cosB / rb / W7 * (-(cosB * y3b / rb + 1) / W7 * W2 + W1 / W7 ^ 2 * W2 * W3 + 1 / 2 *
        W1 / W7 * a / rb ^ 3 * 2 * y3b - a / rb2 + a * y3b / rb2 ^ 2 * 2 * y3b)) /π/(1 - nu)) +
        b2 * (1 / 4 * ((-2 + 2 * nu) * N1 * cotB * ((y3b / rb + 1) / W6 - cosB * W3 / W7) + (2 - 2 * nu) *
        y1 / W6 ^ 2 * W5 * (y3b / rb + 1) + 1 / 2 * (2 - 2 * nu) * y1 / W6 * a / rb ^ 3 * 2 * y3b + (2 -
        2 * nu)* sinB / W7 * W2 - (2 - 2 * nu) * z1b / W7 ^ 2 * W2 * W3 - 1 / 2 * (2 - 2 * nu) * z1b /
        W7 * a / rb ^ 3 * 2 * y3b + 1 / rb * (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb2) - 1 / 2 *
        W8 / rb ^ 3 * (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb2) * 2 * y3b + W8 / rb * (2 * nu *
        y1 / W6 ^ 2 * (y3b / rb + 1) + a * y1 / rb2 ^ 2 * 2 * y3b) -1 / W7 * (cosB * sinB + W1 *
        cotB / rb * ((2 - 2 * nu) * cosB - W1 / W7) + a / rb * (sinB - y3b * z1b / rb2 - z1b *
        W1 / rb / W7)) +W8 / W7 ^ 2 * (cosB * sinB + W1 * cotB / rb * ((2 - 2 * nu) * cosB - W1 /
        W7) +a / rb * (sinB - y3b * z1b / rb2 - z1b * W1 / rb / W7)) * W3 - W8 / W7 * ((cosB *
        y3b / rb + 1)* cotB / rb * ((2 - 2 * nu) * cosB - W1 / W7) - 1 / 2 * W1 * cotB / rb ^ 3 *
        ((2 - 2 * nu) * cosB - W1 / W7) * 2 * y3b + W1 * cotB / rb * (-(cosB * y3b / rb + 1) / W7 +
        W1 / W7 ^ 2 * W3) -1 / 2 * a / rb ^ 3 * (sinB - y3b * z1b / rb2 - z1b * W1 / rb / W7) *
        2 * y3b + a / rb * (-z1b / rb2 - y3b * sinB / rb2 + y3b * z1b / rb2 ^ 2 * 2 * y3b -
        sinB * W1 / rb / W7 - z1b * (cosB * y3b / rb + 1) / rb / W7 + 1 / 2 * z1b * W1 / rb ^ 3 /
        W7 * 2 * y3b + z1b * W1 / rb / W7 ^ 2 * W3))) /π/(1 - nu)) +
        b3 * (1 / 4 * ((2 - 2 * nu) * rFib_ry3 - (2 - 2 * nu) * y2 * sinB / W7 ^ 2 * W2 * W3 - 1 / 2 *
        (2 - 2 * nu) * y2 * sinB / W7 * a / rb ^ 3 * 2 * y3b + y2 * sinB / rb / W7 * (1 + W1 / W7 *
        W2 + a * y3b / rb2) -1 / 2 * y2 * W8 * sinB / rb ^ 3 / W7 * (1 + W1 / W7 * W2 + a * y3b /
        rb2)* 2 * y3b - y2 * W8 * sinB / rb / W7 ^ 2 * (1 + W1 / W7 * W2 + a * y3b / rb2) * W3 +
        y2 * W8 * sinB / rb / W7 * ((cosB * y3b / rb + 1) / W7 * W2 - W1 / W7 ^ 2 * W2 * W3 -
        1 / 2 * W1 / W7 * a / rb ^ 3 * 2 * y3b + a / rb2 - a * y3b / rb2 ^ 2 * 2 * y3b)) /π/(1 - nu)))

    v12 = (b1 / 2 * (1 / 4 * ((-2 + 2 * nu) * N1 * rFib_ry2 * cotB ^ 2 + N1 / W6 * ((1 - W5) * cotB - y1 /
        W6 * W4) -N1 * y2 ^ 2 / W6 ^ 2 * ((1 - W5) * cotB - y1 / W6 * W4) / rb + N1 * y2 / W6 *
        (a / rb ^ 3 * y2 * cotB + y1 / W6 ^ 2 * W4 / rb * y2 + y2 / W6 * a / rb ^ 3 * y1) + N1 *
        cosB * cotB / W7 * W2 - N1 * y2 ^ 2 * cosB * cotB / W7 ^ 2 * W2 / rb - N1 * y2 ^ 2 * cosB *
        cotB / W7 * a / rb ^ 3 + a * W8 * cotB / rb ^ 3 - 3 * a * y2 ^ 2 * W8 * cotB / rb ^ 5 + W8 /
        rb / W6 * (-N1 * cotB + y1 / W6 * W5 + a * y1 / rb2) - y2 ^ 2 * W8 / rb ^ 3 / W6 * (-N1 *
        cotB + y1 / W6 * W5 + a * y1 / rb2) -y2 ^ 2 * W8 / rb2 / W6 ^ 2 * (-N1 * cotB + y1 /
        W6 * W5 + a * y1 / rb2) +y2 * W8 / rb / W6 * (-y1 / W6 ^ 2 * W5 / rb * y2 - y2 / W6 *
        a / rb ^ 3 * y1 - 2 * a * y1 / rb2 ^ 2 * y2) +W8 / rb / W7 * (cosB / W7 * (W1 * (N1 *
        cosB - a / rb)* cotB + (2 - 2 * nu) * (rb * sinB - y1) * cosB) -a * y3b * cosB * cotB /
        rb2) -y2 ^ 2 * W8 / rb ^ 3 / W7 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + (2 -
        2 * nu) * (rb * sinB - y1) * cosB) -a * y3b * cosB * cotB / rb2) -y2 ^ 2 * W8 / rb2 /
        W7 ^ 2 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + (2 - 2 * nu) * (rb * sinB - y1) *
        cosB) -a * y3b * cosB * cotB / rb2) +y2 * W8 / rb / W7 * (-cosB / W7 ^ 2 * (W1 *
        (N1 * cosB - a / rb) * cotB + (2 - 2 * nu) * (rb * sinB - y1) * cosB) /rb *y2+cosB /
        W7 * (1 / rb * cosB * y2 * (N1 * cosB - a / rb) * cotB + W1 * a / rb ^ 3 * y2 * cotB + (2 - 2 *
        nu) /rb*sinB *y2*cosB)+2*a *y3b*cosB*cotB /rb2 ^ 2 * y2)) /π/(1 - nu)) +
        b2 / 2 * (1 / 4 * (N1 * (((2 - 2 * nu) * cotB ^ 2 + nu) / rb * y2 / W6 - ((2 - 2 * nu) * cotB ^ 2 + 1) *
        cosB / rb * y2 / W7) -N1 / W6 ^ 2 * (-N1 * y1 * cotB + nu * y3b - a + a * y1 * cotB / rb +
        y1 ^ 2 / W6 * W4) /rb *y2+N1 /W6 * (-a * y1 * cotB / rb ^ 3 * y2 - y1 ^ 2 / W6 ^
        2 * W4 / rb * y2 - y1 ^ 2 / W6 * a / rb ^ 3 * y2) +N1 * cotB / W7 ^ 2 * (z1b * cosB - a *
        (rb * sinB - y1) / rb / cosB) /rb *y2-N1*cotB /W7 * (-a / rb2 * sinB * y2 /
        cosB + a * (rb * sinB - y1) / rb ^ 3 / cosB * y2) +3 * a * y2 * W8 * cotB / rb ^ 5. * y1 -
        W8 / W6 ^ 2 * (2 * nu + 1 / rb * (N1 * y1 * cotB + a) - y1 ^ 2 / rb / W6 * W5 - a * y1 ^ 2 /
        rb ^ 3) /rb *y2+W8 /W6 * (-1 / rb ^ 3 * (N1 * y1 * cotB + a) * y2 + y1 ^ 2 / rb ^
        3 / W6 * W5 * y2 + y1 ^ 2 / rb2 / W6 ^ 2 * W5 * y2 + y1 ^ 2 / rb2 ^ 2 / W6 * a * y2 + 3 *
        a * y1 ^ 2 / rb ^ 5. * y2) -W8 * cotB / W7 ^ 2 * (-cosB * sinB + a * y1 * y3b / rb ^ 3 /
        cosB + (rb * sinB - y1) / rb * ((2 - 2 * nu) * cosB - W1 / W7 * W9)) /rb *y2+W8*cotB /
        W7 * (-3 * a * y1 * y3b / rb ^ 5. / cosB * y2 + 1 / rb2 * sinB * y2 * ((2 - 2 * nu) * cosB -
        W1 / W7 * W9) -(rb * sinB - y1) / rb ^ 3 * ((2 - 2 * nu) * cosB - W1 / W7 * W9) * y2 + (rb *
        sinB - y1) /rb *(-1 /rb * cosB * y2 / W7 * W9 + W1 / W7 ^ 2 * W9 / rb * y2 + W1 / W7 *
        a / rb ^ 3 / cosB * y2))) /π/(1 - nu)) +
        b3 / 2 * (1 / 4 * (N1 * (1 / W6 * (1 + a / rb) - y2 ^ 2 / W6 ^ 2 * (1 + a / rb) / rb - y2 ^ 2 /
        W6 * a / rb ^ 3 - cosB / W7 * W2 + y2 ^ 2 * cosB / W7 ^ 2 * W2 / rb + y2 ^ 2 * cosB / W7 *
        a / rb ^ 3) -W8 / rb * (a / rb2 + 1 / W6) + y2 ^ 2 * W8 / rb ^ 3 * (a / rb2 + 1 / W6) -
        y2 * W8 / rb * (-2 * a / rb2 ^ 2 * y2 - 1 / W6 ^ 2 / rb * y2) + W8 * cosB / rb / W7 *
        (W1 / W7 * W2 + a * y3b / rb2) - y2 ^ 2 * W8 * cosB / rb ^ 3 / W7 * (W1 / W7 * W2 + a *
        y3b / rb2) -y2 ^ 2 * W8 * cosB / rb2 / W7 ^ 2 * (W1 / W7 * W2 + a * y3b / rb2) + y2 *
        W8 * cosB / rb / W7 * (1 / rb * cosB * y2 / W7 * W2 - W1 / W7 ^ 2 * W2 / rb * y2 - W1 /
        W7 * a / rb ^ 3 * y2 - 2 * a * y3b / rb2 ^ 2 * y2)) /π/(1 - nu)) +
        b1 / 2 * (1 / 4 * (N1 * (((2 - 2 * nu) * cotB ^ 2 - nu) / rb * y1 / W6 - ((2 - 2 * nu) * cotB ^ 2 + 1 -
        2 * nu)* cosB * (y1 / rb - sinB) / W7) +N1 / W6 ^ 2 * (y1 * cotB * (1 - W5) + nu * y3b -
        a + y2 ^ 2 / W6 * W4) /rb *y1-N1 /W6 * ((1 - W5) * cotB + a * y1 ^ 2 * cotB / rb ^ 3 -
        y2 ^ 2 / W6 ^ 2 * W4 / rb * y1 - y2 ^ 2 / W6 * a / rb ^ 3 * y1) -N1 * cosB * cotB / W7 *
        W2 + N1 * z1b * cotB / W7 ^ 2 * W2 * (y1 / rb - sinB) + N1 * z1b * cotB / W7 * a / rb ^
        3 * y1 - a * W8 * cotB / rb ^ 3 + 3 * a * y1 ^ 2 * W8 * cotB / rb ^ 5 - W8 / W6 ^ 2 * (-2 *
        nu + 1 / rb * (N1 * y1 * cotB - a) + y2 ^ 2 / rb / W6 * W5 + a * y2 ^ 2 / rb ^ 3) /rb *
        y1 + W8 / W6 * (-1 / rb ^ 3 * (N1 * y1 * cotB - a) * y1 + 1 / rb * N1 * cotB - y2 ^ 2 /
        rb ^ 3 / W6 * W5 * y1 - y2 ^ 2 / rb2 / W6 ^ 2 * W5 * y1 - y2 ^ 2 / rb2 ^ 2 / W6 * a * y1 -
        3 * a * y2 ^ 2 / rb ^ 5. * y1) -W8 / W7 ^ 2 * (cosB ^ 2 - 1 / rb * (N1 * z1b * cotB + a *
        cosB) +a * y3b * z1b * cotB / rb ^ 3 - 1 / rb / W7 * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB /
        rb * W1)) * (y1 / rb - sinB) + W8 / W7 * (1 / rb ^ 3 * (N1 * z1b * cotB + a * cosB) *
        y1 - 1 / rb * N1 * cosB * cotB + a * y3b * cosB * cotB / rb ^ 3 - 3 * a * y3b * z1b * cotB /
        rb ^ 5. * y1 + 1 / rb ^ 3 / W7 * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB / rb * W1) * y1 + 1 /
        rb / W7 ^ 2 * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB / rb * W1) * (y1 / rb - sinB) - 1 / rb /
        W7 * (-a * cosB * cotB / rb * W1 + a * z1b * cotB / rb ^ 3 * W1 * y1 - a * z1b * cotB /
        rb2 * cosB * y1))) /π/(1 - nu)) +
        b2 / 2 * (1 / 4 * ((2 - 2 * nu) * N1 * rFib_ry1 * cotB ^ 2 - N1 * y2 / W6 ^ 2 * ((W5 - 1) * cotB +
        y1 / W6 * W4) /rb *y1+N1 *y2 /W6 * (-a / rb ^ 3 * y1 * cotB + 1 / W6 * W4 - y1 ^
        2 / W6 ^ 2 * W4 / rb - y1 ^ 2 / W6 * a / rb ^ 3) +N1 * y2 * cotB / W7 ^ 2 * W9 * (y1 /
        rb - sinB) +N1 * y2 * cotB / W7 * a / rb ^ 3 / cosB * y1 + 3 * a * y2 * W8 * cotB / rb ^
        5. * y1 - y2 * W8 / rb ^ 3 / W6 * (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb * (1 / rb + 1 /
        W6)) * y1 - y2 * W8 / rb2 / W6 ^ 2 * (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb * (1 /
        rb + 1 / W6)) * y1 + y2 * W8 / rb / W6 * (-2 * nu / W6 + 2 * nu * y1 ^ 2 / W6 ^ 2 / rb - a /
        rb * (1 / rb + 1 / W6) + a * y1 ^ 2 / rb ^ 3 * (1 / rb + 1 / W6) - a * y1 / rb * (-1 /
        rb ^ 3 * y1 - 1 / W6 ^ 2 / rb * y1)) -y2 * W8 * cotB / rb ^ 3 / W7 * ((-2 + 2 * nu) *
        cosB + W1 / W7 * W9 + a * y3b / rb2 / cosB) * y1 - y2 * W8 * cotB / rb / W7 ^ 2 * ((-2 +
        2 * nu)* cosB + W1 / W7 * W9 + a * y3b / rb2 / cosB) * (y1 / rb - sinB) + y2 * W8 *
        cotB / rb / W7 * (1 / rb * cosB * y1 / W7 * W9 - W1 / W7 ^ 2 * W9 * (y1 / rb - sinB) -
        W1 / W7 * a / rb ^ 3 / cosB * y1 - 2 * a * y3b / rb2 ^ 2 / cosB * y1)) /π/(1 - nu)) +
        b3 / 2 * (1 / 4 * (N1 * (-sinB * (y1 / rb - sinB) / W7 - 1 / W6 * (1 + a / rb) + y1 ^ 2 / W6 ^
        2 * (1 + a / rb) / rb + y1 ^ 2 / W6 * a / rb ^ 3 + cosB / W7 * W2 - z1b / W7 ^ 2 * W2 *
        (y1 / rb - sinB) - z1b / W7 * a / rb ^ 3 * y1) +W8 / rb * (a / rb2 + 1 / W6) - y1 ^ 2 *
        W8 / rb ^ 3 * (a / rb2 + 1 / W6) + y1 * W8 / rb * (-2 * a / rb2 ^ 2 * y1 - 1 / W6 ^ 2 /
        rb * y1) +W8 / W7 ^ 2 * (sinB * (cosB - a / rb) + z1b / rb * (1 + a * y3b / rb2) - 1 /
        rb / W7 * (y2 ^ 2 * cosB * sinB - a * z1b / rb * W1)) * (y1 / rb - sinB) - W8 / W7 *
        (sinB * a / rb ^ 3 * y1 + cosB / rb * (1 + a * y3b / rb2) - z1b / rb ^ 3 * (1 + a * y3b /
        rb2) * y1 - 2 * z1b / rb ^ 5 * a * y3b * y1 + 1 / rb ^ 3 / W7 * (y2 ^ 2 * cosB * sinB - a *
        z1b / rb * W1) * y1 + 1 / rb / W7 ^ 2 * (y2 ^ 2 * cosB * sinB - a * z1b / rb * W1) *
        (y1 / rb - sinB) - 1 / rb / W7 * (-a * cosB / rb * W1 + a * z1b / rb ^ 3 * W1 * y1 - a *
        z1b / rb2 * cosB * y1))) /π/(1 - nu)))

    v13 = (b1 / 2 * (1 / 4 * ((-2 + 2 * nu) * N1 * rFib_ry3 * cotB ^ 2 - N1 * y2 / W6 ^ 2 * ((1 - W5) *
        cotB - y1 / W6 * W4) * (y3b / rb + 1) + N1 * y2 / W6 * (1 / 2 * a / rb ^ 3 * 2 * y3b * cotB +
        y1 / W6 ^ 2 * W4 * (y3b / rb + 1) + 1 / 2 * y1 / W6 * a / rb ^ 3 * 2 * y3b) -N1 * y2 * cosB *
        cotB / W7 ^ 2 * W2 * W3 - 1 / 2 * N1 * y2 * cosB * cotB / W7 * a / rb ^ 3 * 2 * y3b + a /
        rb ^ 3 * y2 * cotB - 3 / 2 * a * y2 * W8 * cotB / rb ^ 5 * 2 * y3b + y2 / rb / W6 * (-N1 *
        cotB + y1 / W6 * W5 + a * y1 / rb2) -1 / 2 * y2 * W8 / rb ^ 3 / W6 * (-N1 * cotB + y1 /
        W6 * W5 + a * y1 / rb2)* 2 * y3b - y2 * W8 / rb / W6 ^ 2 * (-N1 * cotB + y1 / W6 * W5 +
        a * y1 / rb2) * (y3b / rb + 1) + y2 * W8 / rb / W6 * (-y1 / W6 ^ 2 * W5 * (y3b / rb +
        1) -1 / 2 * y1 / W6 * a / rb ^ 3 * 2 * y3b - a * y1 / rb2 ^ 2 * 2 * y3b) +y2 / rb / W7 *
        (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB + (2 - 2 * nu) * (rb * sinB - y1) * cosB) -
        a * y3b * cosB * cotB / rb2) -1 / 2 * y2 * W8 / rb ^ 3 / W7 * (cosB / W7 * (W1 * (N1 *
        cosB - a / rb)* cotB + (2 - 2 * nu) * (rb * sinB - y1) * cosB) -a * y3b * cosB * cotB /
        rb2)* 2 * y3b - y2 * W8 / rb / W7 ^ 2 * (cosB / W7 * (W1 * (N1 * cosB - a / rb) * cotB +
        (2 - 2 * nu) * (rb * sinB - y1) * cosB) -a * y3b * cosB * cotB / rb2) * W3 + y2 * W8 / rb /
        W7 * (-cosB / W7 ^ 2 * (W1 * (N1 * cosB - a / rb) * cotB + (2 - 2 * nu) * (rb * sinB - y1) *
        cosB) * W3 + cosB / W7 * ((cosB * y3b / rb + 1) * (N1 * cosB - a / rb) * cotB + 1 / 2 * W1 *
        a / rb ^ 3 * 2 * y3b * cotB + 1 / 2 * (2 - 2 * nu) / rb * sinB * 2 * y3b * cosB) -a * cosB *
        cotB / rb2 + a * y3b * cosB * cotB / rb2 ^ 2 * 2 * y3b)) /π/(1 - nu)) +
        b2 / 2 * (1 / 4 * (N1 * (((2 - 2 * nu) * cotB ^ 2 + nu) * (y3b / rb + 1) / W6 - ((2 - 2 * nu) * cotB ^
        2 + 1)* cosB * W3 / W7) -N1 / W6 ^ 2 * (-N1 * y1 * cotB + nu * y3b - a + a * y1 * cotB /
        rb + y1 ^ 2 / W6 * W4) * (y3b / rb + 1) + N1 / W6 * (nu - 1 / 2 * a * y1 * cotB / rb ^ 3 * 2 *
        y3b - y1 ^ 2 / W6 ^ 2 * W4 * (y3b / rb + 1) - 1 / 2 * y1 ^ 2 / W6 * a / rb ^ 3 * 2 * y3b) +
        N1 * cotB / W7 ^ 2 * (z1b * cosB - a * (rb * sinB - y1) / rb / cosB) * W3 - N1 * cotB /
        W7 * (cosB * sinB - 1 / 2 * a / rb2 * sinB * 2 * y3b / cosB + 1 / 2 * a * (rb * sinB - y1) /
        rb ^ 3 / cosB * 2 * y3b) -a / rb ^ 3 * y1 * cotB + 3 / 2 * a * y1 * W8 * cotB / rb ^ 5 * 2 *
        y3b + 1 / W6 * (2 * nu + 1 / rb * (N1 * y1 * cotB + a) - y1 ^ 2 / rb / W6 * W5 - a * y1 ^ 2 /
        rb ^ 3) -W8 / W6 ^ 2 * (2 * nu + 1 / rb * (N1 * y1 * cotB + a) - y1 ^ 2 / rb / W6 * W5 - a *
        y1 ^ 2 / rb ^ 3) * (y3b / rb + 1) + W8 / W6 * (-1 / 2 / rb ^ 3 * (N1 * y1 * cotB + a) * 2 *
        y3b + 1 / 2 * y1 ^ 2 / rb ^ 3 / W6 * W5 * 2 * y3b + y1 ^ 2 / rb / W6 ^ 2 * W5 * (y3b / rb +
        1) +1 / 2 * y1 ^ 2 / rb2 ^ 2 / W6 * a * 2 * y3b + 3 / 2 * a * y1 ^ 2 / rb ^ 5 * 2 * y3b) +
        cotB / W7 * (-cosB * sinB + a * y1 * y3b / rb ^ 3 / cosB + (rb * sinB - y1) / rb * ((2 -
        2 * nu)* cosB - W1 / W7 * W9)) -W8 * cotB / W7 ^ 2 * (-cosB * sinB + a * y1 * y3b / rb ^
        3 / cosB + (rb * sinB - y1) / rb * ((2 - 2 * nu) * cosB - W1 / W7 * W9)) * W3 + W8 * cotB /
        W7 * (a / rb ^ 3 / cosB * y1 - 3 / 2 * a * y1 * y3b / rb ^ 5. / cosB * 2 * y3b + 1 / 2 /
        rb2 * sinB * 2 * y3b * ((2 - 2 * nu) * cosB - W1 / W7 * W9) - 1 / 2 * (rb * sinB - y1) / rb ^
        3 * ((2 - 2 * nu) * cosB - W1 / W7 * W9) * 2 * y3b + (rb * sinB - y1) / rb * (-(cosB * y3b /
        rb + 1) /W7 *W9+W1 /W7 ^ 2 * W9 * W3 + 1 / 2 * W1 / W7 * a / rb ^ 3 / cosB * 2 *
        y3b))) /π/(1 - nu)) +
        b3 / 2 * (1 / 4 * (N1 * (-y2 / W6 ^ 2 * (1 + a / rb) * (y3b / rb + 1) - 1 / 2 * y2 / W6 * a /
        rb ^ 3 * 2 * y3b + y2 * cosB / W7 ^ 2 * W2 * W3 + 1 / 2 * y2 * cosB / W7 * a / rb ^ 3 * 2 *
        y3b) -y2 / rb * (a / rb2 + 1 / W6) + 1 / 2 * y2 * W8 / rb ^ 3 * (a / rb2 + 1 / W6) * 2 *
        y3b - y2 * W8 / rb * (-a / rb2 ^ 2 * 2 * y3b - 1 / W6 ^ 2 * (y3b / rb + 1)) + y2 * cosB /
        rb / W7 * (W1 / W7 * W2 + a * y3b / rb2) - 1 / 2 * y2 * W8 * cosB / rb ^ 3 / W7 * (W1 /
        W7 * W2 + a * y3b / rb2)* 2 * y3b - y2 * W8 * cosB / rb / W7 ^ 2 * (W1 / W7 * W2 + a *
        y3b / rb2) * W3 + y2 * W8 * cosB / rb / W7 * ((cosB * y3b / rb + 1) / W7 * W2 - W1 /
        W7 ^ 2 * W2 * W3 - 1 / 2 * W1 / W7 * a / rb ^ 3 * 2 * y3b + a / rb2 - a * y3b / rb2 ^ 2 * 2 *
        y3b)) /π/(1 - nu)) +
        b1 / 2 * (1 / 4 * ((2 - 2 * nu) * (N1 * rFib_ry1 * cotB - y1 / W6 ^ 2 * W5 / rb * y2 - y2 / W6 *
        a / rb ^ 3 * y1 + y2 * cosB / W7 ^ 2 * W2 * (y1 / rb - sinB) + y2 * cosB / W7 * a / rb ^
        3 * y1) -y2 * W8 / rb ^ 3 * (2 * nu / W6 + a / rb2) * y1 + y2 * W8 / rb * (-2 * nu / W6 ^
        2 / rb * y1 - 2 * a / rb2 ^ 2 * y1) -y2 * W8 * cosB / rb ^ 3 / W7 * (1 - 2 * nu - W1 / W7 *
        W2 - a * y3b / rb2) * y1 - y2 * W8 * cosB / rb / W7 ^ 2 * (1 - 2 * nu - W1 / W7 * W2 - a *
        y3b / rb2) * (y1 / rb - sinB) + y2 * W8 * cosB / rb / W7 * (-1 / rb * cosB * y1 / W7 *
        W2 + W1 / W7 ^ 2 * W2 * (y1 / rb - sinB) + W1 / W7 * a / rb ^ 3 * y1 + 2 * a * y3b / rb2 ^
        2 * y1)) /π/(1 - nu)) +
        b2 / 2 * (1 / 4 * ((-2 + 2 * nu) * N1 * cotB * (1 / rb * y1 / W6 - cosB * (y1 / rb - sinB) / W7) -
        (2 - 2 * nu) / W6 * W5 + (2 - 2 * nu) * y1 ^ 2 / W6 ^ 2 * W5 / rb + (2 - 2 * nu) * y1 ^ 2 / W6 *
        a / rb ^ 3 + (2 - 2 * nu) * cosB / W7 * W2 - (2 - 2 * nu) * z1b / W7 ^ 2 * W2 * (y1 / rb -
        sinB) -(2 - 2 * nu) * z1b / W7 * a / rb ^ 3 * y1 - W8 / rb ^ 3 * (N1 * cotB - 2 * nu * y1 /
        W6 - a * y1 / rb2) * y1 + W8 / rb * (-2 * nu / W6 + 2 * nu * y1 ^ 2 / W6 ^ 2 / rb - a / rb2 +
        2 * a * y1 ^ 2 / rb2 ^ 2) +W8 / W7 ^ 2 * (cosB * sinB + W1 * cotB / rb * ((2 - 2 * nu) *
        cosB - W1 / W7) +a / rb * (sinB - y3b * z1b / rb2 - z1b * W1 / rb / W7)) * (y1 / rb -
        sinB) -W8 / W7 * (1 / rb2 * cosB * y1 * cotB * ((2 - 2 * nu) * cosB - W1 / W7) - W1 *
        cotB / rb ^ 3 * ((2 - 2 * nu) * cosB - W1 / W7) * y1 + W1 * cotB / rb * (-1 / rb * cosB *
        y1 / W7 + W1 / W7 ^ 2 * (y1 / rb - sinB)) -a / rb ^ 3 * (sinB - y3b * z1b / rb2 -
        z1b * W1 / rb / W7) * y1 + a / rb * (-y3b * cosB / rb2 + 2 * y3b * z1b / rb2 ^ 2 * y1 -
        cosB * W1 / rb / W7 - z1b / rb2 * cosB * y1 / W7 + z1b * W1 / rb ^ 3 / W7 * y1 + z1b *
        W1 / rb / W7 ^ 2 * (y1 / rb - sinB)))) /π/(1 - nu)) +
        b3 / 2 * (1 / 4 * ((2 - 2 * nu) * rFib_ry1 - (2 - 2 * nu) * y2 * sinB / W7 ^ 2 * W2 * (y1 / rb -
        sinB) -(2 - 2 * nu) * y2 * sinB / W7 * a / rb ^ 3 * y1 - y2 * W8 * sinB / rb ^ 3 / W7 * (1 +
        W1 / W7 * W2 + a * y3b / rb2) * y1 - y2 * W8 * sinB / rb / W7 ^ 2 * (1 + W1 / W7 * W2 +
        a * y3b / rb2) * (y1 / rb - sinB) + y2 * W8 * sinB / rb / W7 * (1 / rb * cosB * y1 /
        W7 * W2 - W1 / W7 ^ 2 * W2 * (y1 / rb - sinB) - W1 / W7 * a / rb ^ 3 * y1 - 2 * a * y3b /
        rb2 ^ 2 * y1)) /π/(1 - nu)))

    v23 = (b1 / 2 * (1 / 4 * (N1 * (((2 - 2 * nu) * cotB ^ 2 - nu) * (y3b / rb + 1) / W6 - ((2 - 2 * nu) *
        cotB ^ 2 + 1 - 2 * nu)* cosB * W3 / W7) +N1 / W6 ^ 2 * (y1 * cotB * (1 - W5) + nu * y3b - a +
        y2 ^ 2 / W6 * W4) * (y3b / rb + 1) - N1 / W6 * (1 / 2 * a * y1 * cotB / rb ^ 3 * 2 * y3b +
        nu - y2 ^ 2 / W6 ^ 2 * W4 * (y3b / rb + 1) - 1 / 2 * y2 ^ 2 / W6 * a / rb ^ 3 * 2 * y3b) -N1 *
        sinB * cotB / W7 * W2 + N1 * z1b * cotB / W7 ^ 2 * W2 * W3 + 1 / 2 * N1 * z1b * cotB / W7 *
        a / rb ^ 3 * 2 * y3b - a / rb ^ 3 * y1 * cotB + 3 / 2 * a * y1 * W8 * cotB / rb ^ 5 * 2 * y3b +
        1 / W6 * (-2 * nu + 1 / rb * (N1 * y1 * cotB - a) + y2 ^ 2 / rb / W6 * W5 + a * y2 ^ 2 /
        rb ^ 3) -W8 / W6 ^ 2 * (-2 * nu + 1 / rb * (N1 * y1 * cotB - a) + y2 ^ 2 / rb / W6 * W5 +
        a * y2 ^ 2 / rb ^ 3) * (y3b / rb + 1) + W8 / W6 * (-1 / 2 / rb ^ 3 * (N1 * y1 * cotB - a) *
        2 * y3b - 1 / 2 * y2 ^ 2 / rb ^ 3 / W6 * W5 * 2 * y3b - y2 ^ 2 / rb / W6 ^ 2 * W5 * (y3b /
        rb + 1) -1 / 2 * y2 ^ 2 / rb2 ^ 2 / W6 * a * 2 * y3b - 3 / 2 * a * y2 ^ 2 / rb ^ 5 * 2 * y3b) +
        1 / W7 * (cosB ^ 2 - 1 / rb * (N1 * z1b * cotB + a * cosB) + a * y3b * z1b * cotB / rb ^
        3 - 1 / rb / W7 * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB / rb * W1)) -W8 / W7 ^ 2 * (cosB ^ 2 -
        1 / rb * (N1 * z1b * cotB + a * cosB) + a * y3b * z1b * cotB / rb ^ 3 - 1 / rb / W7 *
        (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB / rb * W1)) * W3 + W8 / W7 * (1 / 2 / rb ^ 3 * (N1 *
        z1b * cotB + a * cosB)* 2 * y3b - 1 / rb * N1 * sinB * cotB + a * z1b * cotB / rb ^ 3 + a *
        y3b * sinB * cotB / rb ^ 3 - 3 / 2 * a * y3b * z1b * cotB / rb ^ 5 * 2 * y3b + 1 / 2 / rb ^
        3 / W7 * (y2 ^ 2 * cosB ^ 2 - a * z1b * cotB / rb * W1) * 2 * y3b + 1 / rb / W7 ^ 2 * (y2 ^
        2 * cosB ^ 2 - a * z1b * cotB / rb * W1) * W3 - 1 / rb / W7 * (-a * sinB * cotB / rb * W1 +
        1 / 2 * a * z1b * cotB / rb ^ 3 * W1 * 2 * y3b - a * z1b * cotB / rb * (cosB * y3b / rb +
        1)))) /π/(1 - nu)) +
        b2 / 2 * (1 / 4 * ((2 - 2 * nu) * N1 * rFib_ry3 * cotB ^ 2 - N1 * y2 / W6 ^ 2 * ((W5 - 1) * cotB +
        y1 / W6 * W4) * (y3b / rb + 1) + N1 * y2 / W6 * (-1 / 2 * a / rb ^ 3 * 2 * y3b * cotB - y1 /
        W6 ^ 2 * W4 * (y3b / rb + 1) - 1 / 2 * y1 / W6 * a / rb ^ 3 * 2 * y3b) +N1 * y2 * cotB /
        W7 ^ 2 * W9 * W3 + 1 / 2 * N1 * y2 * cotB / W7 * a / rb ^ 3 / cosB * 2 * y3b - a / rb ^ 3 *
        y2 * cotB + 3 / 2 * a * y2 * W8 * cotB / rb ^ 5 * 2 * y3b + y2 / rb / W6 * (N1 * cotB - 2 *
        nu * y1 / W6 - a * y1 / rb * (1 / rb + 1 / W6)) -1 / 2 * y2 * W8 / rb ^ 3 / W6 * (N1 *
        cotB - 2 * nu * y1 / W6 - a * y1 / rb * (1 / rb + 1 / W6))* 2 * y3b - y2 * W8 / rb / W6 ^
        2 * (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb * (1 / rb + 1 / W6)) * (y3b / rb + 1) + y2 *
        W8 / rb / W6 * (2 * nu * y1 / W6 ^ 2 * (y3b / rb + 1) + 1 / 2 * a * y1 / rb ^ 3 * (1 / rb +
        1 / W6)* 2 * y3b - a * y1 / rb * (-1 / 2 / rb ^ 3 * 2 * y3b - 1 / W6 ^ 2 * (y3b / rb +
        1))) +y2 * cotB / rb / W7 * ((-2 + 2 * nu) * cosB + W1 / W7 * W9 + a * y3b / rb2 / cosB) -
        1 / 2 * y2 * W8 * cotB / rb ^ 3 / W7 * ((-2 + 2 * nu) * cosB + W1 / W7 * W9 + a * y3b /
        rb2 / cosB)* 2 * y3b - y2 * W8 * cotB / rb / W7 ^ 2 * ((-2 + 2 * nu) * cosB + W1 / W7 *
        W9 + a * y3b / rb2 / cosB) * W3 + y2 * W8 * cotB / rb / W7 * ((cosB * y3b / rb + 1) /
        W7 * W9 - W1 / W7 ^ 2 * W9 * W3 - 1 / 2 * W1 / W7 * a / rb ^ 3 / cosB * 2 * y3b + a / rb2 /
        cosB - a * y3b / rb2 ^ 2 / cosB * 2 * y3b)) /π/(1 - nu)) +
        b3 / 2 * (1 / 4 * (N1 * (-sinB * W3 / W7 + y1 / W6 ^ 2 * (1 + a / rb) * (y3b / rb + 1) +
        1 / 2 * y1 / W6 * a / rb ^ 3 * 2 * y3b + sinB / W7 * W2 - z1b / W7 ^ 2 * W2 * W3 - 1 / 2 *
        z1b / W7 * a / rb ^ 3 * 2 * y3b) +y1 / rb * (a / rb2 + 1 / W6) - 1 / 2 * y1 * W8 / rb ^
        3 * (a / rb2 + 1 / W6) * 2 * y3b + y1 * W8 / rb * (-a / rb2 ^ 2 * 2 * y3b - 1 / W6 ^ 2 *
        (y3b / rb + 1)) -1 / W7 * (sinB * (cosB - a / rb) + z1b / rb * (1 + a * y3b / rb2) - 1 /
        rb / W7 * (y2 ^ 2 * cosB * sinB - a * z1b / rb * W1)) +W8 / W7 ^ 2 * (sinB * (cosB -
        a / rb) +z1b / rb * (1 + a * y3b / rb2) - 1 / rb / W7 * (y2 ^ 2 * cosB * sinB - a * z1b /
        rb * W1)) * W3 - W8 / W7 * (1 / 2 * sinB * a / rb ^ 3 * 2 * y3b + sinB / rb * (1 + a * y3b /
        rb2) -1 / 2 * z1b / rb ^ 3 * (1 + a * y3b / rb2) * 2 * y3b + z1b / rb * (a / rb2 - a *
        y3b / rb2 ^ 2 * 2 * y3b) +1 / 2 / rb ^ 3 / W7 * (y2 ^ 2 * cosB * sinB - a * z1b / rb *
        W1)* 2 * y3b + 1 / rb / W7 ^ 2 * (y2 ^ 2 * cosB * sinB - a * z1b / rb * W1) * W3 - 1 /
        rb / W7 * (-a * sinB / rb * W1 + 1 / 2 * a * z1b / rb ^ 3 * W1 * 2 * y3b - a * z1b / rb *
        (cosB * y3b / rb + 1)))) /π/(1 - nu)) +
        b1 / 2 * (1 / 4 * ((2 - 2 * nu) * (N1 * rFib_ry2 * cotB + 1 / W6 * W5 - y2 ^ 2 / W6 ^ 2 * W5 /
        rb - y2 ^ 2 / W6 * a / rb ^ 3 - cosB / W7 * W2 + y2 ^ 2 * cosB / W7 ^ 2 * W2 / rb + y2 ^ 2 *
        cosB / W7 * a / rb ^ 3) +W8 / rb * (2 * nu / W6 + a / rb2) - y2 ^ 2 * W8 / rb ^ 3 * (2 *
        nu / W6 + a / rb2) +y2 * W8 / rb * (-2 * nu / W6 ^ 2 / rb * y2 - 2 * a / rb2 ^ 2 * y2) +
        W8 * cosB / rb / W7 * (1 - 2 * nu - W1 / W7 * W2 - a * y3b / rb2) - y2 ^ 2 * W8 * cosB /
        rb ^ 3 / W7 * (1 - 2 * nu - W1 / W7 * W2 - a * y3b / rb2) - y2 ^ 2 * W8 * cosB / rb2 / W7 ^
        2 * (1 - 2 * nu - W1 / W7 * W2 - a * y3b / rb2) + y2 * W8 * cosB / rb / W7 * (-1 / rb *
        cosB * y2 / W7 * W2 + W1 / W7 ^ 2 * W2 / rb * y2 + W1 / W7 * a / rb ^ 3 * y2 + 2 * a *
        y3b / rb2 ^ 2 * y2)) /π/(1 - nu)) +
        b2 / 2 * (1 / 4 * ((-2 + 2 * nu) * N1 * cotB * (1 / rb * y2 / W6 - cosB / rb * y2 / W7) + (2 -
        2 * nu) * y1 / W6 ^ 2 * W5 / rb * y2 + (2 - 2 * nu) * y1 / W6 * a / rb ^ 3 * y2 - (2 - 2 *
        nu) * z1b / W7 ^ 2 * W2 / rb * y2 - (2 - 2 * nu) * z1b / W7 * a / rb ^ 3 * y2 - W8 / rb ^
        3 * (N1 * cotB - 2 * nu * y1 / W6 - a * y1 / rb2) * y2 + W8 / rb * (2 * nu * y1 / W6 ^ 2 /
        rb * y2 + 2 * a * y1 / rb2 ^ 2 * y2) +W8 / W7 ^ 2 * (cosB * sinB + W1 * cotB / rb * ((2 -
        2 * nu)* cosB - W1 / W7) +a / rb * (sinB - y3b * z1b / rb2 - z1b * W1 / rb / W7)) /
        rb * y2 - W8 / W7 * (1 / rb2 * cosB * y2 * cotB * ((2 - 2 * nu) * cosB - W1 / W7) - W1 *
        cotB / rb ^ 3 * ((2 - 2 * nu) * cosB - W1 / W7) * y2 + W1 * cotB / rb * (-cosB / rb *
        y2 / W7 + W1 / W7 ^ 2 / rb * y2) -a / rb ^ 3 * (sinB - y3b * z1b / rb2 - z1b * W1 /
        rb / W7) * y2 + a / rb * (2 * y3b * z1b / rb2 ^ 2 * y2 - z1b / rb2 * cosB * y2 / W7 +
        z1b * W1 / rb ^ 3 / W7 * y2 + z1b * W1 / rb2 / W7 ^ 2 * y2))) /π/(1 - nu)) +
        b3 / 2 * (1 / 4 * ((2 - 2 * nu) * rFib_ry2 + (2 - 2 * nu) * sinB / W7 * W2 - (2 - 2 * nu) * y2 ^ 2 *
        sinB / W7 ^ 2 * W2 / rb - (2 - 2 * nu) * y2 ^ 2 * sinB / W7 * a / rb ^ 3 + W8 * sinB / rb /
        W7 * (1 + W1 / W7 * W2 + a * y3b / rb2) - y2 ^ 2 * W8 * sinB / rb ^ 3 / W7 * (1 + W1 /
        W7 * W2 + a * y3b / rb2) -y2 ^ 2 * W8 * sinB / rb2 / W7 ^ 2 * (1 + W1 / W7 * W2 + a *
        y3b / rb2) +y2 * W8 * sinB / rb / W7 * (1 / rb * cosB * y2 / W7 * W2 - W1 / W7 ^ 2 *
        W2 / rb * y2 - W1 / W7 * a / rb ^ 3 * y2 - 2 * a * y3b / rb2 ^ 2 * y2)) /π/(1 - nu)))

    return v11, v22, v33, v12, v13, v23
end

function strain_tri3_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T) where {T, V}
    A = transform_matrix(P1, P2, P3)
    _strain_tri3_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, λ, μ, A')
end

"""
    stress_tri3_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T) where {T, V}

Compute stress risen from triangular dislocation in elastic *fullspace*.
    Please see [original version (in supporting information)](https://academic.oup.com/gji/article/201/2/1119/572006#86405752)
    for details, especially the **coordinate system** used here.

## Arguments
The same as [`stress_tri3_hs`](@ref)
"""
function stress_tri3_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T) where {T, V}
    exx, eyy, ezz, exy, exz, eyz = strain_tri3_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, λ, μ)
    return _strain2stress(exx, eyy, ezz, exy, exz, eyz, λ, μ)
end

@inline function _strain2stress(exx::T, eyy::T, ezz::T, exy::T, exz::T, eyz::T, λ::T, μ::T) where T
    ekk = exx + eyy + ezz
    σxx = λ * ekk + 2μ * exx
    σyy = λ * ekk + 2μ * eyy
    σzz = λ * ekk + 2μ * ezz
    σxy = 2μ * exy
    σxz = 2μ * exz
    σyz = 2μ * eyz
    return σxx, σyy, σzz, σxy, σxz, σyz
end

@inline function _strain_tri3_fs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T, At::M) where {T, V, M}
    nu = λ / (μ + λ) / 2
    bx, by, bz = Ts, Ss, Ds
    p1, p2, p3 = [zeros(T, 3) for _ in 1: 3]
    x, y, z = coord_trans(X - P2[1], Y - P2[2], Z - P2[3], At)
    p1[1], p1[2], p1[3] = coord_trans(P1[1] - P2[1], P1[2] - P2[2], P1[3] - P2[3], At)
    p3[1], p3[2], p3[3] = coord_trans(P3[1] - P2[1], P3[2] - P2[2], P3[3] - P2[3], At)

    e12 = normalize(p2 - p1)
    e13 = normalize(p3 - p1)
    e23 = normalize(p3 - p2)

    A = acos(dot(e12, e13))
    B = acos(dot(-e12, e23))
    C = acos(dot(e23, e13))

    Trimode = trimodefinder(y, z, x, p1, p2, p3)

    if Trimode == 1
        Exx1Tp, Eyy1Tp, Ezz1Tp, Exy1Tp, Exz1Tp, Eyz1Tp = TDSetupS(x, y, z, A, bx, by, bz, nu, p1, -e13)
        Exx2Tp, Eyy2Tp, Ezz2Tp, Exy2Tp, Exz2Tp, Eyz2Tp = TDSetupS(x, y, z, B, bx, by, bz, nu, p2, e12)
        Exx3Tp, Eyy3Tp, Ezz3Tp, Exy3Tp, Exz3Tp, Eyz3Tp = TDSetupS(x, y, z, C, bx, by, bz, nu, p3, e23)
        exx = Exx1Tp + Exx2Tp + Exx3Tp
        eyy = Eyy1Tp + Eyy2Tp + Eyy3Tp
        ezz = Ezz1Tp + Ezz2Tp + Ezz3Tp
        exy = Exy1Tp + Exy2Tp + Exy3Tp
        exz = Exz1Tp + Exz2Tp + Exz3Tp
        eyz = Eyz1Tp + Eyz2Tp + Eyz3Tp
    end

    if Trimode == -1
        Exx1Tn, Eyy1Tn, Ezz1Tn, Exy1Tn, Exz1Tn, Eyz1Tn = TDSetupS(x, y, z, A, bx, by, bz, nu, p1, e13)
        Exx2Tn, Eyy2Tn, Ezz2Tn, Exy2Tn, Exz2Tn, Eyz2Tn = TDSetupS(x, y, z, B, bx, by, bz, nu, p2, -e12)
        Exx3Tn, Eyy3Tn, Ezz3Tn, Exy3Tn, Exz3Tn, Eyz3Tn = TDSetupS(x, y, z, C, bx, by, bz, nu, p3, -e23)
        exx = Exx1Tn + Exx2Tn + Exx3Tn
        eyy = Eyy1Tn + Eyy2Tn + Eyy3Tn
        ezz = Ezz1Tn + Ezz2Tn + Ezz3Tn
        exy = Exy1Tn + Exy2Tn + Exy3Tn
        exz = Exz1Tn + Exz2Tn + Exz3Tn
        eyz = Eyz1Tn + Eyz2Tn + Eyz3Tn
    end

    if Trimode == 0
        exx = NaN
        eyy = NaN
        ezz = NaN
        exy = NaN
        exz = NaN
        eyz = NaN
    end

    Exx, Eyy, Ezz, Exy, Exz, Eyz = TensTrans(exx, eyy, ezz, exy, exz, eyz, At')
    return Exx, Eyy, Ezz, Exy, Exz, Eyz
end

@inline function TDSetupS(x::T, y::T, z::T, alpha::T, bx::T, by::T, bz::T, nu::T, TriVertex::V, SideVec::V) where {T, V}
    y1 = SideVec[3] * (y - TriVertex[2]) - SideVec[2] * (z - TriVertex[3])
    z1 = SideVec[2] * (y - TriVertex[2]) + SideVec[3] * (z - TriVertex[3])
    by1 = SideVec[3] * by - SideVec[2] * bz
    bz1 = SideVec[2] * by + SideVec[3] * bz

    exx, eyy, ezz, exy, exz, eyz = AngDisStrain(x, y1, z1, -π + alpha, bx, by1, bz1, nu)
    B = zeros(T, 3, 3)
    B[1,1] = one(T)
    B[2,2] = SideVec[3]
    B[2,3] = SideVec[2]
    B[3,2] = -SideVec[2]
    B[3,3] = SideVec[3]
    exx, eyy, ezz, exy, exz, eyz = TensTrans(exx, eyy, ezz, exy, exz, eyz, B)
    return exx, eyy, ezz, exy, exz, eyz
end

function strain_tri3_hs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T) where {T, V}
    A = transform_matrix(P1, P2, P3)
    strms11, strms22, strms33, strms12, strms13, strms23 = _strain_tri3_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, λ, μ, A')
    strfsc11, strfsc22, strfsc33, strfsc12, strfsc13, strfsc23 = _TDstrain_HarFunc(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, λ, μ, A)
    P1[3] *= -one(T)
    P2[3] *= -one(T)
    P3[3] *= -one(T)
    allatsurface = P1[3] == zero(T) && P2[3] == zero(T) && P3[3] == zero(T)
    if !allatsurface
        A *= -one(T)
        A[3,1] *= -one(T)
        A[3,2] *= -one(T)
        A[3,3] *= -one(T)
    end
    stris11, stris22, stris33, stris12, stris13, stris23 = _strain_tri3_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, λ, μ, A')
    if allatsurface
        stris13 *= -one(T)
        stris23 *= -one(T)
    end
    P1[3] *= -one(T)
    P2[3] *= -one(T)
    P3[3] *= -one(T)
    e11 = strms11 + strfsc11 + stris11
    e22 = strms22 + strfsc22 + stris22
    e33 = strms33 + strfsc33 + stris33
    e12 = strms12 + strfsc12 + stris12
    e13 = strms13 + strfsc13 + stris13
    e23 = strms23 + strfsc23 + stris23
    return e11, e22, e33, e12, e13, e23
end

"""
    stress_tri3_hs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T) where {T, V}

Compute stress risen from triangular dislocation in elastic *halfspace*.
    Please see [original version (in supporting information)](https://academic.oup.com/gji/article/201/2/1119/572006#86405752)
    for details, especially the **coordinate system** used here.

## Arguments
- `X`, `Y`, `Z`: observational coordinates
- `P1`, `P2`, `P3`: three triangular vertices coordinates respectively
- `Ss`, `Ds`, `Ts`: triangular dislocation vector, Strike-slip, Dip-slip, Tensile-slip respectively
- `λ`: Lamé's first parameter
- `μ`: shear modulus

## Output
By order: ``σ_{xx}``, ``σ_{yy}``, ``σ_{zz}``,
    ``σ_{xy}``, ``σ_{xz}``, ``σ_{yz}``.
"""
function stress_tri3_hs(X::T, Y::T, Z::T, P1::V, P2::V, P3::V, Ss::T, Ds::T, Ts::T, λ::T, μ::T) where {T, V}
    exx, eyy, ezz, exy, exz, eyz = strain_tri3_hs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, λ, μ)
    return _strain2stress(exx, eyy, ezz, exy, exz, eyz, λ, μ)
end
