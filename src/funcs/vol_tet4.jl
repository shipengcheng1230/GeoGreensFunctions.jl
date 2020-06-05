# Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - May 19, 2018, Los Angeles.
# Translated by Pengcheng Shi (shipengcheng1230@gmail.com) 08/2019

# greens function ϵ vs u for tet4 element
export _disp_vol_tet4, _disp_vol_tet4!
export _strain_vol_tet4, _strain_vol_tet4!
export _stress_vol_tet4, _stress_vol_tet4!

export disp_vol_tet4, disp_vol_tet4!
export stress_vol_tet4, stress_vol_tet4!

function _swap_coord!(x::AbstractVector)
    x[1], x[2], x[3] = x[2], x[1], -x[3]
end

function disp_vol_tet4!(u::AbstractVector{R}, quadrature::Q, x::R, y::R, z::R, A::U, B::U, C::U, D::U, ϵxx::R, ϵxy::R, ϵxz::R, ϵyy::R, ϵyz::R, ϵzz::R, ν::R) where {R, Q, U}
    foreach(x -> _swap_coord!(x), [A, B, C, D])
    _disp_vol_tet4!(u, quadrature, y, x, -z, A, B, C, D, ϵyy, ϵxy, -ϵyz, ϵxx, -ϵxz, ϵzz, ν)
    foreach(x -> _swap_coord!(x), [A, B, C, D])
    u[1], u[2], u[3] = u[2], u[1], -u[3]
end

function disp_vol_tet4(quadrature::Q, x::R, y::R, z::R, A::U, B::U, C::U, D::U, ϵxx::R, ϵxy::R, ϵxz::R, ϵyy::R, ϵyz::R, ϵzz::R, ν::R) where {R, Q, U}
    u = Vector{R}(undef, 3)
    disp_vol_tet4!(u, quadrature, x, y, z, A, B, C, D, ϵxx, ϵxy, ϵxz, ϵyy, ϵyz, ϵzz, ν)
    return u
end

function stress_vol_tet4!(σ::AbstractVector{R}, quadrature::Q, x::R, y::R, z::R, A::U, B::U, C::U, D::U, ϵxx::R, ϵxy::R, ϵxz::R, ϵyy::R, ϵyz::R, ϵzz::R, G::R, ν::R) where {R, Q, U}
    foreach(x -> _swap_coord!(x), [A, B, C, D])
    _stress_vol_tet4!(σ, quadrature, y, x, -z, A, B, C, D, ϵyy, ϵxy, -ϵyz, ϵxx, -ϵxz, ϵzz, G, ν)
    foreach(x -> _swap_coord!(x), [A, B, C, D])
    σ[1], σ[4] = σ[4], σ[1]
    σ[3], σ[5] = -σ[5], -σ[3]
end

function stress_vol_tet4(quadrature::Q, x::R, y::R, z::R, A::U, B::U, C::U, D::U, ϵxx::R, ϵxy::R, ϵxz::R, ϵyy::R, ϵyz::R, ϵzz::R, G::R, ν::R) where {R, Q, U}
    σ = Vector{R}(undef, 6)
    stress_vol_tet4!(σ, quadrature, x, y, z, A, B, C, D, ϵxx, ϵxy, ϵxz, ϵyy, ϵyz, ϵzz, G, ν)
    return σ
end

@doc raw"""
    _disp_vol_tet4(quadrature::Q,
        x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
        e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, nu::R
        ) where {R, U, Q}

Compute displacement arisen from inelastic strain in Tet4 elements.
    Please see [original version](https://bitbucket.org/sbarbot/bssa-2018058/src/default/)
    for complete details, especially the **coordinate system** used here.

## Arguments
- `quadrature`: quadrature rule for integration, see
    [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
- `x1`, `x2`, `x3`: observational position, where ``x_{3} ≥ 0``
- `A`, `B`, `C`, `D`: a list of 3 numbers for each, each of which represents
    coordinates of the vertex. All depth coordinates must be greater or
    equal to 0 (no checking is performed here)
- `e**`: strain components, each is ``ϵ_{11}``, ``ϵ_{12}``, ``ϵ_{13}``,
    ``ϵ_{22}``, ``ϵ_{23}``, ``ϵ_{33}``
- `nu`: poisson ratio

## Output
A vector of 3 numbers, each represents ``u_{1}``, ``u_{2}``, ``u_{3}``

## Notice
- Inplace version: `_disp_vol_tet4!(u, args...)` where u is a vector of 3 numbers.
"""
function _disp_vol_tet4(quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, nu::R
    ) where {R, U, Q}

    u = Vector{R}(undef, 3)
    _disp_vol_tet4!(u, quadrature, x1, x2, x3, A, B, C, D, e11, e12, e13, e22, e23, e33, nu)
    return u
end

function _disp_vol_tet4!(u::W, quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, nu::R
    ) where {R, U, Q, W}

    @assert x3 ≥ zero(R) "Depth must be greater than or equal to 0."
    lambda = 2 * nu / (1 - 2 * nu)
    ekk = e11 + e22 + e33

    nA = cross(C - B, D - B)
    nB = cross(D - C, A - C)
    nC = cross(A - D, B - D)
    nD = cross(B - A, C - A)

    nA /= norm(nA)
    nB /= norm(nB)
    nC /= norm(nC)
    nD /= norm(nD)

    if nA' * (A .- (B .+ C .+ D) ./ 3) > zero(R); nA *= -one(R) end
    if nB' * (B .- (A .+ C .+ D) ./ 3) > zero(R); nB *= -one(R) end
    if nC' * (C .- (B .+ A .+ D) ./ 3) > zero(R); nC *= -one(R) end
    if nD' * (D .- (B .+ C .+ A) ./ 3) > zero(R); nD *= -one(R) end

    ABC = norm(cross(C .- A, B .- A)) / 2
    BCD = norm(cross(D .- B, C .- B)) / 2
    CDA = norm(cross(A .- C, D .- C)) / 2
    DAB = norm(cross(B .- D, A .- D)) / 2

    m11 = lambda * ekk + 2 * e11
    m12 = 2 * e12
    m13 = 2 * e13
    m22 = lambda * ekk + 2 * e22
    m23 = 2 * e23
    m33 = lambda * ekk + 2 * e33

    let lambda=lambda, x1=x1, x2=x2, x3=x3, A=A, B=B, C=C, D=D, e11=e11, e12=e12, e13=e13, e22=e22, e23=e23, e33=e33, nu=nu, nA=nA, nB=nB, nC=nC, ABC=ABC, BCD=BCD, CDA=CDA, DAB=DAB, m11=m11, m12=m12, m13=m13, m22=m22, m23=m23, m33=m33

        function r1(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 - y3)
        end

        function r2(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 + y3)
        end

        function G11(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            (3 - 4 * nu) / r1(y1, y2, y3) + 1 / r2(y1, y2, y3) + (x1 - y1) ^ 2 / r1(y1, y2, y3) ^ 3
            + (3 - 4 * nu) * (x1 - y1) ^ 2 / r2(y1, y2, y3) ^ 3 + 2 * x3 * y3 * (r2(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 5
            + 4 * (1 - 2 * nu) * (1 - nu) * (r2(y1, y2, y3) ^ 2 - (x1 - y1) ^ 2 + r2(y1, y2, y3) * (x3 + y3)) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G12(y1::R, y2::R, y3::R) where R
            ((x1 - y1) * (x2 - y2) / (16 * π * (1 - nu)) * (
            1  / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) / r2(y1, y2, y3) ^ 3 - 6 * x3 * y3 / r2(y1, y2, y3) ^ 5
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G13(y1::R, y2::R, y3::R) where R
            ((x1 - y1) /(16*π*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            -6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 + 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G21(y1::R, y2::R, y3::R) where R
            ((x1 - y1) * (x2 - y2) / (16 * π * (1 - nu)) * (
            1  / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) / r2(y1, y2, y3) ^ 3 - 6 * x3 * y3 / r2(y1, y2, y3) ^ 5
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G22(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            (3 - 4 * nu) / r1(y1, y2, y3) + 1  / r2(y1, y2, y3) + (x2 - y2) ^ 2  / r1(y1, y2, y3) ^ 3
            +(3 - 4 * nu) * (x2 - y2) ^ 2  / r2(y1, y2, y3) ^ 3 + 2 * x3 * y3 * (r2(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 5
            +4 * (1 - 2 * nu) * (1 - nu) * (r2(y1, y2, y3) ^ 2 - (x2 - y2) ^ 2 + r2(y1, y2, y3) * (x3 + y3)) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)))
        end

        function G23(y1::R, y2::R, y3::R) where R
            ((x2 - y2) /(16*π*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            -6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 + 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G31(y1::R, y2::R, y3::R) where R
            ((x1 - y1) /(16*π*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            +6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 - 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G32(y1::R, y2::R, y3::R) where R
            ((x2 - y2) /(16*π*(1-nu)) *(
              (x3 - y3) / r1(y1, y2, y3) ^ 3 + (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 3
            +6 * x3 * y3 * (x3 + y3) / r2(y1, y2, y3) ^ 5 - 4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))))
        end

        function G33(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            (3 - 4 * nu) / r1(y1, y2, y3) + (5 - 12 * nu + 8 * nu ^ 2) / r2(y1, y2, y3) + (x3 - y3) ^ 2  / r1(y1, y2, y3) ^ 3
            +6 * x3 * y3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 5 + ((3 - 4 * nu) * (x3 + y3) ^ 2 - 2 * x3 * y3) / r2(y1, y2, y3) ^ 3))
        end

        function y(u::R, v::R, A::R, B::R, C::R) where R
            A * (1 - u) * (1 - v) / 4 + B * (1 + u) * (1 - v) / 4 + C * (1 + v) / 2
        end

        function IU1(u::R, v::R) where R
            (ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G11(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G21(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G31(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G11(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G21(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G31(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G11(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G21(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G31(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G11(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G21(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G31(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU2(u::R, v::R) where R
            (ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G12(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G22(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G32(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G12(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G22(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G32(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G12(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G22(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G32(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G12(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G22(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G32(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU3(u::R, v::R) where R
            (ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G13(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G23(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G33(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
                + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G13(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G23(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G33(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
                + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G13(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G23(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G33(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
                + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G13(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G23(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
                + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G33(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        fill!(u, zero(R))
        x, w = quadrature
        @inbounds @fastmath for k = 1: length(x)
            @simd for j = 1: length(x)
                u[1] += w[j] * w[k] * (1 - x[k]) * IU1(x[j], x[k])
                u[2] += w[j] * w[k] * (1 - x[k]) * IU2(x[j], x[k])
                u[3] += w[j] * w[k] * (1 - x[k]) * IU3(x[j], x[k])
            end
        end

        return nothing
    end
end

@doc raw"""
    _stress_vol_tet4(quadrature::Q,
        x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
        e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, G::R, nu::R
        ) where {R, U, Q}

Compute stress arisen from inelastic strain in Tet4 elements.
    Please see [original version](https://bitbucket.org/sbarbot/bssa-2018058/src/default/)
    for complete details, especially the **coordinate system** used here.

## Arguments
- `quadrature`: quadrature rule for integration, see
    [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
- `x1`, `x2`, `x3`: observational position
- `A`, `B`, `C`, `D`: a list of 3 numbers for each, each of which represents
    coordinates of the vertex
- `e**`: strain components, each is ``ϵ_{11}``, ``ϵ_{12}``, ``ϵ_{13}``,
    ``ϵ_{22}``, ``ϵ_{23}``, ``ϵ_{33}``
- `G`: shear modulus
- `nu`: poisson ratio

## Output
A vector of 6 numbers, each represents ``σ_{11}``, ``σ_{12}``, ``σ_{13}``,
    ``σ_{22}``, ``σ_{23}``, ``σ_{33}``

## Notice
- Inplace version: `_stress_vol_tet4!(u, args...)` where u is a vector of 6 numbers.
"""
function _stress_vol_tet4(quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, G::R, nu::R
    ) where {R, U, Q}

    u = Vector{R}(undef, 6)
    _stress_vol_tet4!(u, quadrature, x1, x2, x3, A, B, C, D, e11, e12, e13, e22, e23, e33, G, nu)
    return u
end

function _stress_vol_tet4!(u::W, quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, G::R, nu::R
    ) where {R, U, Q, W}

    _strain_vol_tet4!(u, quadrature, x1, x2, x3, A, B, C, D, e11, e12, e13, e22, e23, e33, G, nu)

    lambda = 2 * nu / (1 - 2 * nu)
    ekk = u[1] + u[4] + u[6]
    u[1] = 2G * u[1] + G * lambda * ekk
    u[2] *= 2G
    u[3] *= 2G
    u[4] = 2G * u[4] + G * lambda * ekk
    u[5] *= 2G
    u[6] = 2G * u[6] + G * lambda * ekk

    return nothing
end

function _strain_vol_tet4(quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, G::R, nu::R
    ) where {R, U, Q}

    u = Vector{R}(undef, 6)
    _strain_vol_tet4!(u, quadrature, x1, x2, x3, A, B, C, D, e11, e12, e13, e22, e23, e33, G, nu)
    return u
end

function _strain_vol_tet4!(u::W, quadrature::Q,
    x1::R, x2::R, x3::R, A::U, B::U, C::U, D::U,
    e11::R, e12::R, e13::R, e22::R, e23::R, e33::R, G::R, nu::R
    ) where {R, U, Q, W}

    @assert x3 ≥ zero(R) "Depth must be greater than or equal to 0."
    lambda = 2 * nu / (1 - 2 * nu)
    ekk = e11 + e22 + e33

    nA = cross(C - B, D - B)
    nB = cross(D - C, A - C)
    nC = cross(A - D, B - D)
    nD = cross(B - A, C - A)

    nA /= norm(nA)
    nB /= norm(nB)
    nC /= norm(nC)
    nD /= norm(nD)

    if nA' * (A .- (B .+ C .+ D) ./ 3) > zero(R); nA *= -one(R) end
    if nB' * (B .- (A .+ C .+ D) ./ 3) > zero(R); nB *= -one(R) end
    if nC' * (C .- (B .+ A .+ D) ./ 3) > zero(R); nC *= -one(R) end
    if nD' * (D .- (B .+ C .+ A) ./ 3) > zero(R); nD *= -one(R) end

    ABC = norm(cross(C .- A, B .- A)) / 2
    BCD = norm(cross(D .- B, C .- B)) / 2
    CDA = norm(cross(A .- C, D .- C)) / 2
    DAB = norm(cross(B .- D, A .- D)) / 2

    m11 = lambda * ekk + 2 * e11
    m12 = 2 * e12
    m13 = 2 * e13
    m22 = lambda * ekk + 2 * e22
    m23 = 2 * e23
    m33 = lambda * ekk + 2 * e33

    let lambda=lambda, x1=x1, x2=x2, x3=x3, A=A, B=B, C=C, D=D, e11=e11, e12=e12, e13=e13, e22=e22, e23=e23, e33=e33, G=G, nu=nu, nA=nA, nB=nB, nC=nC, ABC=ABC, BCD=BCD, CDA=CDA, DAB=DAB, m11=m11, m12=m12, m13=m13, m22=m22, m23=m23, m33=m33

        function r1(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 - y3)
        end

        function r2(y1::R, y2::R, y3::R) where R
            hypot(x1 - y1, x2 - y2, x3 + y3)
        end

        function y(u::R, v::R, A::R, B::R, C::R) where R
            A * (1 - u) * (1 - v) / 4 + B * (1 + u) * (1 - v) / 4 + C * (1 + v) / 2
        end

        function G11d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (
            -(3 - 4 * nu) / r1(y1, y2, y3) ^ 3
            -1  / r2(y1, y2, y3) ^ 3
            +(2 * r1(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (2 * r2(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (3 * r2(y1, y2, y3) ^ 2 - 5 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            -8 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - 2 * nu) * (1 - nu) * (x1 - y1) ^ 2  / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +8 * (1 - 2 * nu) * (1 - nu) * (x1 - y1) ^ 2  / (r2(y1, y2, y3) ^ 2  * (r2(y1, y2, y3) + x3 + y3) ^ 3) );)
        end

        function G11d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x2 - y2) * (
            -(3 - 4 * nu) / r1(y1, y2, y3) ^ 3
            -1  / r2(y1, y2, y3) ^ 3
            -3 * (x1 - y1) ^ 2  / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x1 - y1) ^ 2  / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (r2(y1, y2, y3) ^ 2 - 5 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - 2 * nu) * (1 - nu) * (x1 - y1) ^ 2  * (3 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 3) ))
        end

        function G11d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            -(3 - 4 * nu) * (x3 - y3) / r1(y1, y2, y3) ^ 3
            -(x3 + y3) / r2(y1, y2, y3) ^ 3
            -3 * (x1 - y1) ^ 2  * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x1 - y1) ^ 2  * (x3 + y3) / r2(y1, y2, y3) ^ 5
            +2 * y3 * (r2(y1, y2, y3) ^ 2 - 3 * x3 * (x3 + y3)) / r2(y1, y2, y3) ^ 5
            -6 * y3 * (x1 - y1) ^ 2  * (r2(y1, y2, y3) ^ 2 - 5 * x3 * (x3 + y3)) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))
            +4 * (1 - 2 * nu) * (1 - nu) * (x1 - y1) ^ 2  * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G21d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x2 - y2) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (r2(y1, y2, y3) ^ 2 - 5 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - 2 * nu) * (1 - nu) * (x1 - y1) ^ 2  * (3 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 3) ))
        end

        function G21d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (r2(y1, y2, y3) ^ 2 - 5 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - 2 * nu) * (1 - nu) * (x2 - y2) ^ 2  * (3 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 3) ))
        end

        function G21d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (x2 - y2) * (
            -3 * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x3 + y3) / r2(y1, y2, y3) ^ 5
            -6 * y3 * (r2(y1, y2, y3) ^ 2 - 5 * x3 * (x3 + y3)) / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G31d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            +(x3 - y3) * (r1(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (x3 - y3) * (r2(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 5
            +6 * y3 * x3 * (x3 + y3) * (r2(y1, y2, y3) ^ 2 - 5 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))
            +4 * (1 - 2 * nu) * (1 - nu) * (x1 - y1) ^ 2  * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G31d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (x2 - y2) * (
            -3 * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 5
            -30 * y3 * x3 * (x3 + y3) / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G31d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x3 - y3) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x3 ^ 2 - y3 ^ 2)) / r2(y1, y2, y3) ^ 5
            +6 * y3 * (2 * x3 + y3) / r2(y1, y2, y3) ^ 5
            -30 * y3 * x3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) / r2(y1, y2, y3) ^ 3 ))
        end

        function G12d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x2 - y2) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (r2(y1, y2, y3) ^ 2 - 5 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - nu) * (1 - 2 * nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - nu) * (1 - 2 * nu) * (x1 - y1) ^ 2  * (3 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 3) ))
        end

        function G12d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (r2(y1, y2, y3) ^ 2 - 5 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - nu) * (1 - 2 * nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - nu) * (1 - 2 * nu) * (x2 - y2) ^ 2  * (3 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 3) ))
        end

        function G12d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (x2 - y2) * (
            -3 * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x3 + y3) / r2(y1, y2, y3) ^ 5
            -6 * y3 * (r2(y1, y2, y3) ^ 2 - 5 * x3 * (x3 + y3)) / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G22d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (
            -(3 - 4 * nu) / r1(y1, y2, y3) ^ 3
            -1  / r2(y1, y2, y3) ^ 3
            -3 * (x2 - y2) ^ 2  / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x2 - y2) ^ 2  / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (r2(y1, y2, y3) ^ 2 - 5 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - 2 * nu) * (1 - nu) * (x2 - y2) ^ 2  * (3 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 3) ))
        end

        function G22d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x2 - y2) * (
            -(3 - 4 * nu) / r1(y1, y2, y3) ^ 3
            -1  / r2(y1, y2, y3) ^ 3
            +(2 * r1(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (2 * r2(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (3 * r2(y1, y2, y3) ^ 2 - 5 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 7
            -12 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3) ^ 2)
            +4 * (1 - 2 * nu) * (1 - nu) * (x2 - y2) ^ 2  * (3 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 3) ))
        end

        function G22d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            -(3 - 4 * nu) * (x3 - y3) / r1(y1, y2, y3) ^ 3
            -(x3 + y3) / r2(y1, y2, y3) ^ 3
            -3 * (x2 - y2) ^ 2  * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x2 - y2) ^ 2  * (x3 + y3) / r2(y1, y2, y3) ^ 5
            +2 * y3 * (r2(y1, y2, y3) ^ 2 - 3 * x3 * (x3 + y3)) / r2(y1, y2, y3) ^ 5
            -6 * y3 * (x2 - y2) ^ 2  * (r2(y1, y2, y3) ^ 2 - 5 * x3 * (x3 + y3)) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))
            +4 * (1 - 2 * nu) * (1 - nu) * (x2 - y2) ^ 2  * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G32d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (x2 - y2) * (
            -3 * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 5
            -30 * y3 * x3 * (x3 + y3) / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G32d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            +(x3 - y3) * (r1(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (x3 - y3) * (r2(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 5
            +6 * y3 * x3 * (x3 + y3) * (r2(y1, y2, y3) ^ 2 - 5 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))
            +4 * (1 - 2 * nu) * (1 - nu) * (x2 - y2) ^ 2  * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G32d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x2 - y2) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x3 - y3) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x3 ^ 2 - y3 ^ 2)) / r2(y1, y2, y3) ^ 5
            +6 * y3 * (2 * x3 + y3) / r2(y1, y2, y3) ^ 5
            -30 * y3 * x3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) / r2(y1, y2, y3) ^ 3 ))
        end

        function G13d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            +(x3 - y3) * (r1(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (x3 - y3) * (r2(y1, y2, y3) ^ 2 - 3 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (x3 + y3) * (r2(y1, y2, y3) ^ 2 - 5 * (x1 - y1) ^ 2) / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))
            -4 * (1 - 2 * nu) * (1 - nu) * (x1 - y1) ^ 2  * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G13d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (x2 - y2) * (
            -3 * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 5
            +30 * y3 * x3 * (x3 + y3) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G13d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x3 - y3) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x3 ^ 2 - y3 ^ 2)) / r2(y1, y2, y3) ^ 5
            -6 * y3 * (2 * x3 + y3) / r2(y1, y2, y3) ^ 5
            +30 * y3 * x3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / r2(y1, y2, y3) ^ 3 ))
        end

        function G23d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (x2 - y2) * (
            -3 * (x3 - y3) / r1(y1, y2, y3) ^ 5
            -3 * (3 - 4 * nu) * (x3 - y3) / r2(y1, y2, y3) ^ 5
            +30 * y3 * x3 * (x3 + y3) / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G23d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            +(x3 - y3) * (r1(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (x3 - y3) * (r2(y1, y2, y3) ^ 2 - 3 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 5
            -6 * y3 * x3 * (x3 + y3) * (r2(y1, y2, y3) ^ 2 - 5 * (x2 - y2) ^ 2) / r2(y1, y2, y3) ^ 7
            +4 * (1 - 2 * nu) * (1 - nu) / (r2(y1, y2, y3) * (r2(y1, y2, y3) + x3 + y3))
            -4 * (1 - 2 * nu) * (1 - nu) * (x2 - y2) ^ 2  * (2 * r2(y1, y2, y3) + x3 + y3) / (r2(y1, y2, y3) ^ 3  * (r2(y1, y2, y3) + x3 + y3) ^ 2) ))
        end

        function G23d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x2 - y2) * (
            +(r1(y1, y2, y3) ^ 2 - 3 * (x3 - y3) ^ 2) / r1(y1, y2, y3) ^ 5
            +(3 - 4 * nu) * (r2(y1, y2, y3) ^ 2 - 3 * (x3 ^ 2 - y3 ^ 2)) / r2(y1, y2, y3) ^ 5
            -6 * y3 * (2 * x3 + y3) / r2(y1, y2, y3) ^ 5
            +30 * y3 * x3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 7
            -4 * (1 - 2 * nu) * (1 - nu) / r2(y1, y2, y3) ^ 3 ))
        end

        function G33d1(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x1 - y1) * (
            -(3 - 4 * nu) / r1(y1, y2, y3) ^ 3
            -(5 - 12 * nu + 8 * nu ^ 2) / r2(y1, y2, y3) ^ 3
            -3 * (x3 - y3) ^ 2  / r1(y1, y2, y3) ^ 5
            -30 * y3 * x3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 7
            -3 * (3 - 4 * nu) * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 5
            +6 * y3 * x3 / r2(y1, y2, y3) ^ 5 ))
        end

        function G33d2(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (x2 - y2) * (
            -(3 - 4 * nu) / r1(y1, y2, y3) ^ 3
            -(5 - 12 * nu + 8 * nu ^ 2) / r2(y1, y2, y3) ^ 3
            -3 * (x3 - y3) ^ 2  / r1(y1, y2, y3) ^ 5
            -30 * y3 * x3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 7
            -3 * (3 - 4 * nu) * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 5
            +6 * y3 * x3 / r2(y1, y2, y3) ^ 5 ))
        end

        function G33d3(y1::R, y2::R, y3::R) where R
            (1 / (16 * π * (1 - nu)) * (
            -(3 - 4 * nu) * (x3 - y3) / r1(y1, y2, y3) ^ 3
            -(5 - 12 * nu + 8 * nu ^ 2) * (x3 + y3) / r2(y1, y2, y3) ^ 3
            +(x3 - y3) * (2 * r1(y1, y2, y3) ^ 2 - 3 * (x3 - y3) ^ 2) / r1(y1, y2, y3) ^ 5
            +6 * y3 * (x3 + y3) ^ 2  / r2(y1, y2, y3) ^ 5
            +6 * y3 * x3 * (x3 + y3) * (2 * r2(y1, y2, y3) ^ 2 - 5 * (x3 + y3) ^ 2) / r2(y1, y2, y3) ^ 7
            +(3 - 4 * nu) * (x3 + y3) * (2 * r2(y1, y2, y3) ^ 2 - 3 * (x3 + y3) ^ 2) / r2(y1, y2, y3) ^ 5
            -2 * y3 * (r2(y1, y2, y3) ^ 2 - 3 * x3 * (x3 + y3)) / r2(y1, y2, y3) ^ 5 ))
        end

        function IU11(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G11d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G21d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G31d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G11d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G21d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G31d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G11d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G21d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G31d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G11d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G21d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G31d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU12(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G11d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G21d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G31d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G11d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G21d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G31d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G11d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G21d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G31d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G11d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G21d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G31d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU13(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G11d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G21d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G31d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G11d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G21d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G31d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G11d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G21d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G31d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G11d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G21d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G31d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU21(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G12d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G22d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G32d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G12d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G22d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G32d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G12d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G22d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G32d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G12d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G22d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G32d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU22(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G12d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G22d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G32d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G12d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G22d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G32d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G12d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G22d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G32d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G12d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G22d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G32d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU23(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G12d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G22d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G32d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G12d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G22d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G32d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G12d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G22d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G32d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G12d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G22d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G32d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU31(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G13d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G23d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G33d1(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G13d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G23d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G33d1(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G13d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G23d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G33d1(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G13d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G23d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G33d1(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU32(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G13d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G23d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G33d2(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G13d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G23d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G33d2(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G13d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G23d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G33d2(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G13d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G23d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G33d2(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        function IU33(u::R, v::R) where R
            (+ ABC / 4 * (m11 * nD[1] + m12 * nD[2] + m13 * nD[3]) * G13d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m12 * nD[1] + m22 * nD[2] + m23 * nD[3]) * G23d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + ABC / 4 * (m13 * nD[1] + m23 * nD[2] + m33 * nD[3]) * G33d3(y(u, v, A[1], B[1], C[1]), y(u, v, A[2], B[2], C[2]), y(u, v, A[3], B[3], C[3]))
            + BCD / 4 * (m11 * nA[1] + m12 * nA[2] + m13 * nA[3]) * G13d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m12 * nA[1] + m22 * nA[2] + m23 * nA[3]) * G23d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + BCD / 4 * (m13 * nA[1] + m23 * nA[2] + m33 * nA[3]) * G33d3(y(u, v, B[1], C[1], D[1]), y(u, v, B[2], C[2], D[2]), y(u, v, B[3], C[3], D[3]))
            + CDA / 4 * (m11 * nB[1] + m12 * nB[2] + m13 * nB[3]) * G13d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m12 * nB[1] + m22 * nB[2] + m23 * nB[3]) * G23d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + CDA / 4 * (m13 * nB[1] + m23 * nB[2] + m33 * nB[3]) * G33d3(y(u, v, C[1], D[1], A[1]), y(u, v, C[2], D[2], A[2]), y(u, v, C[3], D[3], A[3]))
            + DAB / 4 * (m11 * nC[1] + m12 * nC[2] + m13 * nC[3]) * G13d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m12 * nC[1] + m22 * nC[2] + m23 * nC[3]) * G23d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3]))
            + DAB / 4 * (m13 * nC[1] + m23 * nC[2] + m33 * nC[3]) * G33d3(y(u, v, D[1], A[1], B[1]), y(u, v, D[2], A[2], B[2]), y(u, v, D[3], A[3], B[3])))
        end

        u11, u12, u13, u21, u22, u23, u31, u32, u33 = zero(R), zero(R), zero(R), zero(R), zero(R), zero(R), zero(R), zero(R), zero(R)

        x, w = quadrature
        @inbounds @fastmath for k = 1: length(x)
            @simd for j = 1: length(x)
                c = w[j] * w[k] * (1 - x[k])
                u11 += c * IU11(x[j], x[k])
                u12 += c * IU12(x[j], x[k])
                u13 += c * IU13(x[j], x[k])
                u21 += c * IU21(x[j], x[k])
                u22 += c * IU22(x[j], x[k])
                u23 += c * IU23(x[j], x[k])
                u31 += c * IU31(x[j], x[k])
                u32 += c * IU32(x[j], x[k])
                u33 += c * IU33(x[j], x[k])
            end
        end

        omega =
            (heaviside(((A[1] + B[1] + C[1]) / 3 - x1) * nD[1] + ((A[2] + B[2] + C[2]) / 3 - x2) * nD[2] + ((A[3] + B[3] + C[3]) / 3 - x3) * nD[3])
            * heaviside(((B[1] + C[1] + D[1]) / 3 - x1) * nA[1] + ((B[2] + C[2] + D[2]) / 3 - x2) * nA[2] + ((B[3] + C[3] + D[3]) / 3 - x3) * nA[3])
            * heaviside(((C[1] + D[1] + A[1]) / 3 - x1) * nB[1] + ((C[2] + D[2] + A[2]) / 3 - x2) * nB[2] + ((C[3] + D[3] + A[3]) / 3 - x3) * nB[3])
            * heaviside(((D[1] + A[1] + B[1]) / 3 - x1) * nC[1] + ((D[2] + A[2] + B[2]) / 3 - x2) * nC[2] + ((D[3] + A[3] + B[3]) / 3 - x3) * nC[3]))

        u[1] = u11 - omega * e11
        u[2] = (u12 + u21) / 2 - omega * e12
        u[3] = (u13 + u31) / 2 - omega * e13
        u[4] = u22 - omega * e22
        u[5] = (u23 + u32) / 2 - omega * e23
        u[6] = u33 - omega * e33

        return nothing
    end
end
