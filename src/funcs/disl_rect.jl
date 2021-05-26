# Reference journal article:
# Nikkhoo, M., Walter, T. R., Lundgren, P. R., Prats-Iraola, P. (2016):
# Compound dislocation models (CDMs) for volcano deformation analyses.
# Submitted to Geophysical Journal International
#
# Copyright (c) 2016 Mehdi Nikkhoo
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
# created: 2013.1.25
# Last modified: 2016.10.18
#
# Section 2.1, Physics of Earthquakes and Volcanoes
# Department 2, Geophysics
# Helmholtz Centre Potsdam
# German Research Centre for Geosciences (GFZ)
#
# email:
# mehdi.nikkhoo@gfz-potsdam.de
# mehdi.nikkhoo@gmail.com
#
# website:
# http://www.volcanodeformation.com
#
# Translated by Pengcheng Shi (shipengcheng1230@gmail.com) 12/2018
# The broadcast here causes tremendous amount of compiliation time, needs rewriting!

# export disp_rect_fs, strain_rect_fs, stress_rect_fs

# using LinearAlgebra

# function disp_rect_fs(x::V, y::V, z::V, x₀::T, y₀::T, depth::T, L::T, W::T, plunge::T, dip::T, strike::T, rake::T, slip::T, opening::T, nu::T, ref::S) where {V<:AbstractVector, T, S}
#     @assert size(x) == size(y) == size(z) "Unmatched length of observational coordinates."

#     bx = opening
#     by = slip * cosd(rake)
#     bz = slip * sind(rake)
#     Rz1 = [cosd(plunge) sind(plunge) zero(T); -sind(plunge) cosd(plunge) zero(T); zero(T) zero(T) one(T)]
#     Ry = [cosd(dip) zero(T) sind(dip); zero(T) one(T) zero(T); -sind(dip) zero(T) cosd(dip)]
#     Rz2 = [cosd(strike) sind(strike) zero(T); -sind(strike) cosd(strike) zero(T); zero(T) zero(T) one(T)]
#     Rt = Rz2 * Ry * Rz1

#     Pt1 = [-W/2, L/2, zero(T)]
#     Pt2 = [-W/2, -L/2, zero(T)]
#     Pt3 = [W/2, -L/2, zero(T)]
#     Pt4 = [W/2, L/2, zero(T)]

#     Ptr = _ref_pts(ref, L, W)
#     Pr = [x₀, y₀, -depth] - Rt * Ptr
#     P1 = Rt * Pt1 + Pr
#     P2 = Rt * Pt2 + Pr
#     P3 = Rt * Pt3 + Pr
#     P4 = Rt * Pt4 + Pr

#     ez = [zero(T), zero(T), one(T)]
#     Vnorm = Rt * ez
#     Vstrike = [sind(strike), cosd(strike), zero(T)]
#     Vdip = Vnorm × Vstrike

#     Pm = (P1 + P2 + P3 + P4) / 4
#     p1, p2, p3, p4 = [zeros(T, 3) for _ in 1: 4]
#     At = hcat(Vnorm, Vstrike, Vdip)'

#     x, y, z = coord_trans(x .- Pm[1], y .- Pm[2], z .- Pm[3], At)
#     p1[1], p1[2], p1[3] = coord_trans(P1[1] - Pm[1], P1[2] - Pm[2], P1[3] - Pm[3], At)
#     p2[1], p2[2], p2[3] = coord_trans(P2[1] - Pm[1], P2[2] - Pm[2], P2[3] - Pm[3], At)
#     p3[1], p3[2], p3[3] = coord_trans(P3[1] - Pm[1], P3[2] - Pm[2], P3[3] - Pm[3], At)
#     p4[1], p4[2], p4[3] = coord_trans(P4[1] - Pm[1], P4[2] - Pm[2], P4[3] - Pm[3], At)

#     e12 = normalize(p2 - p1)
#     e23 = normalize(p3 - p2)
#     e34 = normalize(p4 - p3)
#     e14 = normalize(p4 - p1)

#     Rectmode = @views rectmodefInder(y, z, x, p1[2:3], p2[2:3], p3[2:3], p4[2:3])

#     # use list comprehension causes being unable infer type (until Julia v1.2), see one example mentioned in the link below
#     # https://discourse.julialang.org/t/comprehension-type-inference/4324
#     # u, v, w = [zeros(T, length(x)) for _ ∈ 1: 3]
#     u = zeros(T, length(x)); v = zeros(T, length(x)); w = zeros(T, length(x))

#     casepLog = Rectmode .== 1
#     casenLog = Rectmode .== -1
#     casezLog = Rectmode .== 0

#     if any(casepLog)
#         xp, yp, zp = x[casepLog], y[casepLog], z[casepLog]
#         utmp, vtmp, wtmp = [Vector{T}(undef, length(xp)) for _ ∈ 1: 3]

#         utmp, vtmp, wtmp = RDSetupD(xp, yp, zp, bx, by, bz, nu, p1, e14)
#         u[casepLog] .+= utmp; v[casepLog] .+= vtmp; w[casepLog] .+= wtmp

#         utmp, vtmp, wtmp = RDSetupD(xp, yp, zp, bx, by, bz, nu, p2, -e12)
#         u[casepLog] .+= utmp; v[casepLog] .+= vtmp; w[casepLog] .+= wtmp

#         utmp, vtmp, wtmp = RDSetupD(xp, yp, zp, bx, by, bz, nu, p3, -e23)
#         u[casepLog] .+= utmp; v[casepLog] .+= vtmp; w[casepLog] .+= wtmp

#         utmp, vtmp, wtmp = RDSetupD(xp, yp, zp, bx, by, bz, nu, p4, -e34)
#         u[casepLog] .+= utmp; v[casepLog] .+= vtmp; w[casepLog] .+= wtmp
#     end
#     if any(casenLog)
#         xn, yn, zn = x[casenLog], y[casenLog], z[casenLog]
#         utmp, vtmp, wtmp = [Vector{T}(undef, length(xn)) for _ ∈ 1: 3]

#         utmp, vtmp, wtmp = RDSetupD(xn, yn, zn, bx, by, bz, nu, p1, -e14)
#         u[casenLog] .+= utmp; v[casenLog] .+= vtmp; w[casenLog] .+= wtmp

#         utmp, vtmp, wtmp = RDSetupD(xn, yn, zn, bx, by, bz, nu, p2, e12)
#         u[casenLog] .+= utmp; v[casenLog] .+= vtmp; w[casenLog] .+= wtmp

#         utmp, vtmp, wtmp = RDSetupD(xn, yn, zn, bx, by, bz, nu, p3, e23)
#         u[casenLog] .+= utmp; v[casenLog] .+= vtmp; w[casenLog] .+= wtmp

#         utmp, vtmp, wtmp = RDSetupD(xn, yn, zn, bx, by, bz, nu, p4, e34)
#         u[casenLog] .+= utmp; v[casenLog] .+= vtmp; w[casenLog] .+= wtmp
#     end
#     if any(casezLog)
#         u[casezLog] .= NaN; v[casezLog] .= NaN; w[casezLog] .= NaN
#     end

#     Fi = BurgersFuncRD(x, y, z, p1, p2, p3, p4)
#     @. u += Fi * bx
#     @. v += Fi * by
#     @. w += Fi * bz
#     u, v, w = coord_trans(u, v, w, At')
#     return u, v, w
# end

# function rd_disp_hs()
#     #TODO
# end

# function strain_rect_fs(x::V, y::V, z::V, x₀::T, y₀::T, depth::T, L::T, W::T, plunge::T, dip::T, strike::T, rake::T, slip::T, opening::T, λ::T, μ::T, ref::S) where {V<:AbstractVector, T, S}
#     nu = λ / (λ + μ) / 2

#     # lots of repeated codes here
#     bx = opening
#     by = slip * cosd(rake)
#     bz = slip * sind(rake)
#     Rz1 = [cosd(plunge) sind(plunge) zero(T); -sind(plunge) cosd(plunge) zero(T); zero(T) zero(T) one(T)]
#     Ry = [cosd(dip) zero(T) sind(dip); zero(T) one(T) zero(T); -sind(dip) zero(T) cosd(dip)]
#     Rz2 = [cosd(strike) sind(strike) zero(T); -sind(strike) cosd(strike) zero(T); zero(T) zero(T) one(T)]
#     Rt = Rz2 * Ry * Rz1

#     Pt1 = [-W/2, L/2, zero(T)]
#     Pt2 = [-W/2, -L/2, zero(T)]
#     Pt3 = [W/2, -L/2, zero(T)]
#     Pt4 = [W/2, L/2, zero(T)]

#     Ptr = _ref_pts(ref, L, W)
#     Pr = [x₀, y₀, -depth] - Rt * Ptr
#     P1 = Rt * Pt1 + Pr
#     P2 = Rt * Pt2 + Pr
#     P3 = Rt * Pt3 + Pr
#     P4 = Rt * Pt4 + Pr

#     ez = [zero(T), zero(T), one(T)]
#     Vnorm = Rt * ez
#     Vstrike = [sind(strike), cosd(strike), zero(T)]
#     Vdip = Vnorm × Vstrike

#     Pm = (P1 + P2 + P3 + P4) / 4
#     p1, p2, p3, p4 = [zeros(T, 3) for _ in 1: 4]
#     At = hcat(Vnorm, Vstrike, Vdip)'

#     x, y, z = coord_trans(x .- Pm[1], y .- Pm[2], z .- Pm[3], At)
#     p1[1], p1[2], p1[3] = coord_trans(P1[1] - Pm[1], P1[2] - Pm[2], P1[3] - Pm[3], At)
#     p2[1], p2[2], p2[3] = coord_trans(P2[1] - Pm[1], P2[2] - Pm[2], P2[3] - Pm[3], At)
#     p3[1], p3[2], p3[3] = coord_trans(P3[1] - Pm[1], P3[2] - Pm[2], P3[3] - Pm[3], At)
#     p4[1], p4[2], p4[3] = coord_trans(P4[1] - Pm[1], P4[2] - Pm[2], P4[3] - Pm[3], At)

#     e12 = normalize(p2 - p1)
#     e23 = normalize(p3 - p2)
#     e34 = normalize(p4 - p3)
#     e14 = normalize(p4 - p1)

#     Rectmode = @views rectmodefInder(y, z, x, p1[2:3], p2[2:3], p3[2:3], p4[2:3])

#     casepLog = Rectmode .== 1
#     casenLog = Rectmode .== -1
#     casezLog = Rectmode .== 0
#     lenx = length(x)

#     # use list comprehension causes being unable infer type (until Julia v1.2), see one example mentioned in the link below
#     # https://discourse.julialang.org/t/comprehension-type-inference/4324
#     # exx, eyy, ezz, exy, exz, eyz = [zeros(T, length(x)) for _ in 1: 6]
#     exx, eyy, ezz, exy, exz, eyz = zeros(T, lenx), zeros(T, lenx), zeros(T, lenx), zeros(T, lenx), zeros(T, lenx), zeros(T, lenx)

#     if any(casepLog)
#         xp, yp, zp = x[casepLog], y[casepLog], z[casepLog]
#         exxt, eyyt, ezzt, exyt, exzt, eyzt = [Vector{T}(undef, length(xp)) for _ ∈ 1: 6]

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xp, yp, zp, bx, by, bz, nu, p1, e14)
#         exx[casepLog] .+= exxt; eyy[casepLog] .+= eyyt; ezz[casepLog] .+= ezzt; exy[casepLog] .+= exyt; exz[casepLog] .+= exzt; eyz[casepLog] .+= eyzt

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xp, yp, zp, bx, by, bz, nu, p2, -e12)
#         exx[casepLog] .+= exxt; eyy[casepLog] .+= eyyt; ezz[casepLog] .+= ezzt; exy[casepLog] .+= exyt; exz[casepLog] .+= exzt; eyz[casepLog] .+= eyzt

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xp, yp, zp, bx, by, bz, nu, p3, -e23)
#         exx[casepLog] .+= exxt; eyy[casepLog] .+= eyyt; ezz[casepLog] .+= ezzt; exy[casepLog] .+= exyt; exz[casepLog] .+= exzt; eyz[casepLog] .+= eyzt

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xp, yp, zp, bx, by, bz, nu, p4, -e34)
#         exx[casepLog] .+= exxt; eyy[casepLog] .+= eyyt; ezz[casepLog] .+= ezzt; exy[casepLog] .+= exyt; exz[casepLog] .+= exzt; eyz[casepLog] .+= eyzt
#     end
#     if any(casenLog)
#         xn, yn, zn = x[casenLog], y[casenLog], z[casenLog]
#         exxt, eyyt, ezzt, exyt, exzt, eyzt = [Vector{T}(undef, length(xn)) for _ ∈ 1: 6]

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xn, yn, zn, bx, by, bz, nu, p1, -e14)
#         exx[casenLog] .+= exxt; eyy[casenLog] .+= eyyt; ezz[casenLog] .+= ezzt; exy[casenLog] .+= exyt; exz[casenLog] .+= exzt; eyz[casenLog] .+= eyzt

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xn, yn, zn, bx, by, bz, nu, p2, e12)
#         exx[casenLog] .+= exxt; eyy[casenLog] .+= eyyt; ezz[casenLog] .+= ezzt; exy[casenLog] .+= exyt; exz[casenLog] .+= exzt; eyz[casenLog] .+= eyzt

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xn, yn, zn, bx, by, bz, nu, p3, e23)
#         exx[casenLog] .+= exxt; eyy[casenLog] .+= eyyt; ezz[casenLog] .+= ezzt; exy[casenLog] .+= exyt; exz[casenLog] .+= exzt; eyz[casenLog] .+= eyzt

#         exxt, eyyt, ezzt, exyt, exzt, eyzt = RDSetupS(xn, yn, zn, bx, by, bz, nu, p4, e34)
#         exx[casenLog] .+= exxt; eyy[casenLog] .+= eyyt; ezz[casenLog] .+= ezzt; exy[casenLog] .+= exyt; exz[casenLog] .+= exzt; eyz[casenLog] .+= eyzt
#     end
#     if any(casezLog)
#         exx[casezLog] .= NaN; eyy[casezLog] .= NaN; ezz[casezLog] .= NaN; exy[casezLog] .= NaN; exz[casezLog] .= NaN; eyz[casezLog] .= NaN
#     end
#     exx, eyy, ezz, exy, exz, eyz = TensTrans(exx, eyy, ezz, exy, exz, eyz, At')
# end

# function rd_strain_hs()
#     #TODO
# end

# function stress_rect_fs(x::V, y::V, z::V, x₀::T, y₀::T, depth::T, L::T, W::T, plunge::T, dip::T, strike::T, rake::T, slip::T, opening::T, λ::T, μ::T, ref::S) where {V, T, S}
#     exx, eyy, ezz, exy, exz, eyz = strain_rect_fs(x, y, z, x₀, y₀, depth, L, W, plunge, dip, strike, rake, slip, opening, λ, μ, ref)
#     ekk = exx + eyy + ezz
#     Sxx = @. 2μ * exx + λ * ekk
#     Syy = @. 2μ * eyy + λ * ekk
#     Szz = @. 2μ * ezz + λ * ekk
#     Sxy = @. 2μ * exy
#     Sxz = @. 2μ * exz
#     Syz = @. 2μ * eyz
#     return Sxx, Syy, Szz, Sxy, Sxz, Syz
# end

# function rd_stress_hs()
#     #TODO
# end

# @inline _ref_pts(::Val{:pc}, L::T, W::T) where T = zeros(T, 3)
# @inline _ref_pts(::Val{:p1}, L::T, W::T) where T = [-W/2, L/2, zero(T)]
# @inline _ref_pts(::Val{:p2}, L::T, W::T) where T = [-W/2, -L/2, zero(T)]
# @inline _ref_pts(::Val{:p3}, L::T, W::T) where T = [W/2, -L/2, zero(T)]
# @inline _ref_pts(::Val{:p4}, L::T, W::T) where T = [W/2, L/2, zero(T)]
# @inline _ref_pts(::Val{:mp12}, L::T, W::T) where T = [-W/2, zero(T), zero(T)]
# @inline _ref_pts(::Val{:mp23}, L::T, W::T) where T = [zero(T), -L/2, zero(T)]
# @inline _ref_pts(::Val{:mp34}, L::T, W::T) where T = [W/2, zero(T), zero(T)]
# @inline _ref_pts(::Val{:mp41}, L::T, W::T) where T = [zero(T), L/2, zero(T)]
# @inline _ref_pts(::Val{:mp21}, L::T, W::T) where T = _ref_pts(Val(:mp12), L, W)
# @inline _ref_pts(::Val{:mp32}, L::T, W::T) where T = _ref_pts(Val(:mp23), L, W)
# @inline _ref_pts(::Val{:mp43}, L::T, W::T) where T = _ref_pts(Val(:mp34), L, W)
# @inline _ref_pts(::Val{:mp14}, L::T, W::T) where T = _ref_pts(Val(:mp14), L, W)

# @inline function rectmodefInder(x::T, y::T, z::T, p1::U, p2::U, p3::U, p4::U) where {T, U}
#     pm = (p1 + p2 + p3 + p4) / 4
#     e21 = normalize(p1 - p2)
#     e41 = normalize(p1 - p4)
#     A = hcat(e21, e41)
#     rectmode = ones(Int, length(x))
#     r = A' * hcat(p1 - pm, p2 - pm, p3 - pm, p4 - pm)
#     P1, P2, P3, P4 = @views r[:,1], r[:,2], r[:,3], r[:,4]

#     @inbounds @simd for i ∈ eachindex(rectmode)
#         _x, _y = x[i] - pm[1], y[i] - pm[2]
#         r1 = A[1,1] * _x + A[2,1] * _y
#         r2 = A[1,2] * _x + A[2,2] * _y
#         _x, _y = r1, r2
#         if (_x ≥ 0 && _y ≥ 0 && ((_y - P1[2]) < (_x - P1[1]))) ||
#             (_x ≤ 0 && _y ≥ 0 && ((_y - P2[2]) > -(_x - P2[1]))) ||
#             (_x ≤ 0 && _y ≤ 0 && ((_y - P3[2]) > (_x - P3[1]))) ||
#             (_x ≥ 0 && _y ≤ 0 && ((_y - P4[2]) < -(_x - P4[1]))) ||
#             (_x < P1[1] && _x > P3[1] && _y < P1[2] && _y > P3[2])
#             rectmode[i] = -1
#         elseif ((_x == P1[1] || _x == P3[1]) && _y ≤ P1[2] && _y ≥ P3[2] && z[i] == 0) ||
#             ((_y == P1[2] || _y == P3[2]) && _x ≤ P1[1] && _x ≥ P3[1] && z[i] == 0)
#             rectmode[i] = 0
#         end
#     end
#     rectmode
# end

# @inline function RDSetupD(x::U, y::U, z::U, bx::T, by::T, bz::T, nu::T, RDVertex::V, SideVec::V) where {T, U, V}
#     y1 = @. SideVec[3] * (y - RDVertex[2]) - SideVec[2] * (z - RDVertex[3])
#     z1 = @. SideVec[2] * (y - RDVertex[2]) + SideVec[3] * (z - RDVertex[3])
#     by1 = @. SideVec[3] * by - SideVec[2] * bz
#     bz1 = @. SideVec[2] * by + SideVec[3] * bz

#     u, v0, w0 = AngDisDisp(x, y1, z1, -π/2, bx, by1, bz1, nu)
#     v = @. SideVec[3] * v0 + SideVec[2] * w0
#     w = @. -SideVec[2] * v0 + SideVec[3] * w0
#     return u, v, w
# end

# @inline function BurgersFuncRD(x::T, y::T, z::T, p1::U, p2::U, p3::U, p4::U) where {T, U}
#     Fi = zeros(size(x))
#     Ind = @. (abs(y) ≤ abs(p1[2] - p2[2])/2) & (abs(z) ≤ abs(p1[3] - p4[3])/2)

#     if any(Ind)
#         xI, yI, zI = x[Ind], y[Ind], z[Ind]

#         FiD1 = @. √(xI ^ 2 + (yI - p1[2]) ^ 2 + (zI - p1[3]) ^ 2) - (zI - p1[3]) - (yI - p1[2])
#         FiD2 = @. √(xI ^ 2 + (yI - p2[2]) ^ 2 + (zI - p2[3]) ^ 2) - (zI - p2[3]) - (yI - p2[2])
#         FiD3 = @. √(xI ^ 2 + (yI - p3[2]) ^ 2 + (zI - p3[3]) ^ 2) - (zI - p3[3]) - (yI - p3[2])
#         FiD4 = @. √(xI ^ 2 + (yI - p4[2]) ^ 2 + (zI - p4[3]) ^ 2) - (zI - p4[3]) - (yI - p4[2])

#         FiNt = @. xI * (FiD1 + FiD3) * (FiD2 * FiD4 - xI ^ 2) - xI * (FiD2 + FiD4) * (FiD1 * FiD3 - xI ^ 2)
#         FiDt = @. (FiD1 * FiD3 - xI ^ 2) * (FiD2 * FiD4 - xI ^ 2) + xI ^ 2 * (FiD1 + FiD3) * (FiD2 + FiD4)
#         Fi[Ind] = @. atan(FiNt, FiDt) / 2 / π
#     elseif !all(Ind)
#         bool = .!Ind
#         xO, yO, zO = x[bool], y[bool], z[bool]
#         a = cat(-xO, p1[2] .- yO, p1[3] .- zO; dims=2)
#         b = cat(-xO, p2[2] .- yO, p2[3] .- zO; dims=2)
#         c = cat(-xO, p3[2] .- yO, p3[3] .- zO; dims=2)
#         d = cat(-xO, p4[2] .- yO, p4[3] .- zO; dims=2)
#         na, nb, nc, nd = vec(.√(sum(abs2, a; dims=2))), vec(.√(sum(abs2, b; dims=2))), vec(.√(sum(abs2, c; dims=2))), vec(.√(sum(abs2, d; dims=2)))
#         @views begin
#             FiN1t = @. a[:,1] * (b[:,2] * c[:,3] - b[:,3] * c[:,2]) - a[:,2] * (b[:,1] * c[:,3]-b[:,3] * c[:,1]) + a[:,3] * (b[:,1] * c[:,2]-b[:,2] * c[:,1])
#             FiD1t = @. na * nb * nc + $vec($sum(a * b; dims=2)) * nc + $vec($sum(a * c; dims=2)) * nb + $vec($sum(b * c; dims=2)) * na
#             FiN2t = @. a[:,1] * (c[:,2] * d[:,3] - c[:,3] * d[:,2]) - a[:,2] * (c[:,1] * d[:,3] - c[:,3] * d[:,1]) + a[:,3] * (c[:,1] * d[:,2]-c[:,2] * d[:,1])
#             FiD2t = @. na * nc * nd + $vec($sum(a * c; dims=2)) * nd + $vec($sum(a * d; dims=2)) * nc + $vec($sum(c * d; dims=2)) * na
#         end
#         Fi[bool] = @. -2 * atan(FiN1t * FiD2t + FiN2t * FiD1t, FiD1t * FiD2t - FiN1t * FiN2t) / 4 / π
#     end
#     return Fi
# end

# @inline function RDSetupS(x::T, y::T, z::T, bx::U, by::U, bz::U, nu::U, RDVertex::V, SideVec::V) where {T, U, V}
#     y1 = @. SideVec[3] * (y - RDVertex[2]) - SideVec[2] * (z - RDVertex[3])
#     z1 = @. SideVec[2] * (y - RDVertex[2]) + SideVec[3] * (z - RDVertex[3])
#     by1 = @. SideVec[3] * by - SideVec[2] * bz
#     bz1 = @. SideVec[2] * by + SideVec[3] * bz

#     exx, eyy, ezz, exy, exz, eyz = AngDisStrain(x, y1, z1, -pi/2, bx, by1, bz1, nu)

#     B = [1 0 0; 0 SideVec[3] SideVec[2]; 0 -SideVec[2] SideVec[3]]
#     exx, eyy, ezz, exy, exz, eyz = TensTrans(exx, eyy, ezz, exy, exz, eyz, B)
#     return exx, eyy, ezz, exy, exz, eyz
# end

# function disp_rect_fs(x::V, y::V, z::V, x₀::T, y₀::T, depth::T, L::T, W::T, plunge::T, dip::T, strike::T, rake::T, slip::T, opening::T, nu::T, ref::S) where {V<:Number, T, S}
#     u, v, w = disp_rect_fs([x], [y], [z], x₀, y₀, depth, L, W, plunge, dip, strike, rake, slip, opening, nu, ref)
#     u[1], v[1], w[1]
# end

# function strain_rect_fs(x::V, y::V, z::V, x₀::T, y₀::T, depth::T, L::T, W::T, plunge::T, dip::T, strike::T, rake::T, slip::T, opening::T, λ::T, μ::T, ref::S) where {V<:Number, T, S}
#     out = strain_rect_fs([x], [y], [z], x₀, y₀, depth, L, W, plunge, dip, strike, rake, slip, opening, λ, μ, ref)
#     [q[1] for q in out]
# end
