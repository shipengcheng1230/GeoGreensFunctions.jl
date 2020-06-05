# Author: Yoshimitsu Okada (https://www.researchgate.net/profile/Yoshimitsu_Okada)
# Translated by Pengcheng Shi (shipengcheng1230@gmail.com) 06/2018

# An example wrapper for DC3D in julia as below:
# function dc3d_fortran(x::T, y::T, z::T, α::T, dep::T, dip::T, al1::T, al2::T, aw1::T, aw2::T,
#     disl1::T, disl2::T, disl3::T) where {T <: AbstractFloat}
#
#     # initial return values
#     # `RefValue{T}` may be also viable other than `Array{T, 1}`
#     ux = Array{Float64}(1)
#     uy = Array{Float64}(1)
#     uz = Array{Float64}(1)
#     uxx = Array{Float64}(1)
#     uyx = Array{Float64}(1)
#     uzx = Array{Float64}(1)
#     uxy = Array{Float64}(1)
#     uyy = Array{Float64}(1)
#     uzy = Array{Float64}(1)
#     uxz = Array{Float64}(1)
#     uyz = Array{Float64}(1)
#     uzz = Array{Float64}(1)
#     iret = Array{Int64}(1)
#
#     # call okada's code which is renamed as "__dc3d__" (see binding rename shown below)
#     # input args tuple must be syntactically written instead of a variable assigned
#     # macros could be used to simplify this in the future
#     ccall((:__dc3d__, "dc3d.so"), Void,
#         (
#             Ref{Float64},
#             Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
#             Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
#             Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
#             Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
#             Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
#             Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
#             Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}}, Ptr{Array{Float64,1}},
#             Ref{Int64},
#         ),
#         α, x, y, z, dep, dip, al1, al2, aw1, aw2, disl1, disl2, disl3,
#         ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz,
#         iret,
#     )
#
#     # results valid iff iret[1] == 0
#     return (
#         iret[1],
#         ux[1], uy[1], uz[1],
#         uxx[1], uyx[1], uzx[1],
#         uxy[1], uyy[1], uzy[1],
#         uxz[1], uyz[1], uzz[1]
#     )
# end
# ```
#
#
# The corresponding fortran module is:
# ```fortran
# MODULE okada
#   USE, INTRINSIC :: iso_c_binding
#   IMPLICIT NONE
# CONTAINS
#
#   SUBROUTINE dc3d_wrapper(&
#        & alpha, &
#        & x, y, z, &
#        & depth, dip, &
#        & al1, al2, &
#        & aw1, aw2, &
#        & disl1, disl2, disl3, &
#        & ux, uy, uz, &
#        & uxx, uyx, uzx, &
#        & uxy, uyy, uzy, &
#        & uxz, uyz, uzz, &
#        & iret) BIND(C, NAME='__dc3d__')
#
#     REAL*8 :: &
#          & alpha, &
#          & x, y, z, &
#          & depth, dip, &
#          & al1, al2, &
#          & aw1, aw2, &
#          & disl1, disl2, disl3, &
#          & ux, uy, uz, &
#          & uxx, uyx, uzx, &
#          & uxy, uyy, uzy, &
#          & uxz, uyz, uzz
#
#     INTEGER*8 :: iret
#
#     CALL dc3d(&
#          & alpha, &
#          & x, y, z, &
#          & depth, dip, &
#          & al1, al2, &
#          & aw1, aw2, &
#          & disl1, disl2, disl3, &
#          & ux, uy, uz, &
#          & uxx, uyx, uzx, &
#          & uxy, uyy, uzy, &
#          & uxz, uyz, uzz, &
#          & iret)
#
#   END SUBROUTINE dc3d_wrapper
#
# END MODULE okada
# ```
#
#
# A sample of makefile is as below:
# ```make
# # Build Okada's code for calculating deformation due to a fault model
# #
# CC = gfortran
# CFLAGS = -fPIC -w -O3
# LDFLAGS = -shared
#
# SRCS = dc3d.f okada.f90
# OBJS = \$(SRCS:.c=.o)
#
# TARGET = dc3d.so
#
# \$(TARGET): \$(OBJS)
#     \$(CC) \$(CFLAGS) \$(LDFLAGS) -o \$(TARGET) \$(OBJS)
#

export dc3d

"""
    dc3d(x::T, y::T, z::T, α::T, dep::T, dip::T,
        al::Union{A, SubArray}, aw::Union{A, SubArray}, disl::A
        ) where {T <: Number, A <: AbstractVecOrMat{T}}

Calculate displacements and gradient of displacements due to a rectangular dislocation
    in an elastic isotropic halfspace.

Please see [dc3d](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) for details.
    Also this fault coordinate system is widely used in this package.

## Arguments
- `x`, `y`, `z`: observational position, where ``z ≤ 0``
- `α`: elastic constant
- `dep`: depth of fault origin
- `dip`: dip angle in degree
- `al`: a vector of 2 numbers, indicating along strike (x-axis) spanning
- `aw`: a vector of 2 numbers, indicating along downdip (y-z plane) spanning
- `disl`: a vector of 3 numbers, indicating dislocation in along-strike,
    along-downdip and tensile respectively.

## Output
A vector of 12 numbers, each is ``u_{x}``, ``u_{y}``, ``u_{z}``, ``u_{x,x}``,
    ``u_{y,x}``, ``u_{z,x}``, ``u_{x,y}``, ``u_{y,y}``, ``u_{z,y}``, ``u_{x,z}``
    ``u_{y,z}``, ``u_{z,z}``.
"""
function dc3d(x::T, y::T, z::T, α::T, dep::T, dip::T, al::A, aw::A, disl::AbstractVector{<:Real}) where {T<:Real, A<:AbstractVecOrMat{<:Real}}

    z > zero(T) && return zeros(T, 12)

    u, du, dua, dub, duc = [zeros(T, 12) for _ in 1: 5]

    xi = x .- al
    d = dep + z
    sc1 = shared_constants_1(α, dip)
    sd, cd = sc1[6], sc1[7]
    p = y * cd + d * sd
    q = y * sd - d * cd
    et = p .- aw

    if q ≈ zero(T) && ((xi[1] * xi[2] ≤ zero(T) && et[1] * et[2] ≈ zero(T)) ||	(et[1] * et[2] ≤ zero(T) && xi[1] * xi[2] ≈ zero(T)))
        return zeros(T, 12)
    end

    kxi, ket = [zeros(T, 2) for _ in 1: 2]
    r12 = hypot(xi[1], et[2], q)
    r21 = hypot(xi[2], et[1], q)
    r22 = hypot(xi[2], et[2], q)

    (xi[1] ≤ zero(T) && (r21 + xi[2]) ≈ zero(T)) && (kxi[1] = one(T))
    (xi[1] ≤ zero(T) && (r22 + xi[2]) ≈ zero(T)) && (kxi[2] = one(T))
    (et[1] ≤ zero(T) && (r12 + et[2]) ≈ zero(T)) && (ket[1] = one(T))
    (et[1] ≤ zero(T) && (r22 + et[2]) ≈ zero(T)) && (ket[2] = one(T))

    for k = 1: 2, j = 1: 2
        sc2 = shared_constants_2(xi[j], et[k], q, sd, cd, kxi[k], ket[j])
        dua = ua(sc1, sc2, xi[j], et[k], q, disl)
        for i = 1: 3: 10
            du[i] = -dua[i]
            du[i+1] = -dua[i+1] * cd + dua[i+2] * sd
            du[i+2] = -dua[i+1] * sd - dua[i+2] * cd
            i < 10 && continue
            du[i] = -du[i]
            du[i+1] = -du[i+1]
            du[i+2] = -du[i+2]
        end
        for i = 1: 12
            if j + k ≠ 3  u[i] += du[i] end
            if j + k == 3  u[i] -= du[i] end
        end
    end

    d = dep - z
    p = y * cd + d * sd
    q = y * sd - d * cd
    et = p .- aw

    if q ≈ zero(T) && ((xi[1] * xi[2] ≤ zero(T) && et[1] * et[2] ≈ zero(T)) ||	(et[1] * et[2] ≤ zero(T) && xi[1] * xi[2] ≈ zero(T)))
        return zeros(T, 12)
    end

    kxi, ket = [zeros(T, 2) for _ in 1: 2]
    r12 = hypot(xi[1], et[2], q)
    r21 = hypot(xi[2], et[1], q)
    r22 = hypot(xi[2], et[2], q)

    (xi[1] ≤ zero(T) && (r21 + xi[2]) ≈ zero(T)) && (kxi[1] = one(T))
    (xi[1] ≤ zero(T) && (r22 + xi[2]) ≈ zero(T)) && (kxi[2] = one(T))
    (et[1] ≤ zero(T) && (r12 + et[2]) ≈ zero(T)) && (ket[1] = one(T))
    (et[1] ≤ zero(T) && (r22 + et[2]) ≈ zero(T)) && (ket[2] = one(T))

    for k = 1: 2, j = 1: 2
        sc2 = shared_constants_2(xi[j], et[k], q, sd, cd, kxi[k], ket[j])
        dua = ua(sc1, sc2, xi[j], et[k], q, disl)
        dub = ub(sc1, sc2, xi[j], et[k], q, disl)
        duc = uc(sc1, sc2, xi[j], et[k], q, z,  disl)
        for i = 1: 3: 10
            du[i] = dua[i] + dub[i] + z * duc[i]
            du[i+1] = (dua[i+1] + dub[i+1] + z * duc[i+1]) * cd - (dua[i+2] + dub[i+2] + z * duc[i+2]) * sd
            du[i+2] = (dua[i+1] + dub[i+1] - z * duc[i+1]) * sd + (dua[i+2] + dub[i+2] - z * duc[i+2]) * cd
            i < 10 && continue
            du[10] += duc[1]
            du[11] += duc[2] * cd - duc[3] * sd
            du[12] -= duc[2] * sd + duc[3] * cd
        end
        for i = 1: 12
            j + k ≠ 3 && (u[i] += du[i])
            j + k == 3 && (u[i] -= du[i])
        end
    end
    return u
end

function ua(sc1::B1, sc2::B2, xi::T, et::T, q::T, disl::A
    ) where {T <: Number, A <: AbstractArray{T}, B1 <: NTuple{12, T}, B2 <: NTuple{24, T}}
    alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d = sc1
    xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz = sc2

    u, du = zeros(T, 12), zeros(T, 12)

    xy = xi * y11
    qx = q * x11
    qy = q * y11

    if disl[1] ≉ zero(T)
        du[1] = tt / 2 + alp2 * xi * qy
        du[2] = alp2 * q / r
        du[3] = alp1 * ale - alp2 * q * qy
        du[4] = -alp1 * qy - alp2 * xi2 * q * y32
        du[5] = -alp2 * xi * q / r3
        du[6] = alp1 * xy + alp2 * xi * q2 * y32
        du[7] = alp1 * xy * sd + alp2 * xi * fy + d / 2 * x11
        du[8] = alp2 * ey
        du[9] = alp1 * (cd / r + qy * sd) - alp2 * q * fy
        du[10] = alp1 * xy * cd + alp2 * xi * fz + y / 2 * x11
        du[11] = alp2 * ez
        du[12] = -alp1 * (sd / r - qy * cd) - alp2 * q * fz
        u .+= disl[1] / 2π * du
    end

    if disl[2] ≉ zero(T)
        du[1] = alp2 * q / r
        du[2] = tt / 2 + alp2 * et * qx
        du[3] = alp1 * alx - alp2 * q * qx
        du[4] = -alp2 * xi * q / r3
        du[5] = -qy / 2 - alp2 * et * q / r3
        du[6] = alp1 / r + alp2 * q2 / r3
        du[7] = alp2 * ey
        du[8] = alp1 * d * x11 + xy / 2 * sd + alp2 * et * gy
        du[9] = alp1 * y * x11 - alp2 * q * gy
        du[10] = alp2 * ez
        du[11] = alp1 * y * x11 + xy / 2 * cd + alp2 * et * gz
        du[12] = -alp1 * d * x11 - alp2 * q * gz
        u .+= disl[2] / 2π * du
    end

    if disl[3] ≉ zero(T)
        du[1] = -alp1 * ale - alp2 * q * qy
        du[2] = -alp1 * alx - alp2 * q * qx
        du[3] = tt / 2 - alp2 * (et * qx + xi * qy)
        du[4] = -alp1 * xy + alp2 * xi * q2 * y32
        du[5] = -alp1 / r + alp2 * q2 / r3
        du[6] = -alp1 * qy - alp2 * q * q2 * y32
        du[7] = -alp1 * (cd / r + qy * sd) - alp2 * q * fy
        du[8] = -alp1 * y * x11 - alp2 * q * gy
        du[9] = alp1 * (d * x11 + xy * sd) + alp2 * q * hy
        du[10] = alp1 * (sd / r - qy * cd) - alp2 * q * fz
        du[11] = alp1 * d * x11 - alp2 * q * gz
        du[12] = alp1 * (y * x11 + xy * cd) + alp2 * q * hz
        u .+= disl[3] / 2π * du
    end
    return u
end

function ub(sc1::B1, sc2::B2, xi::T, et::T, q::T, disl::A
    ) where {T <: Number, A <: AbstractArray{T}, B1 <: NTuple{12, T}, B2 <: NTuple{24, T}}
    alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d = sc1
    xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz = sc2

    u, du = zeros(T, 12), zeros(T, 12)

    rd = r + d
    d11 = one(T) / (r * rd)
    aj2 = xi * y / rd * d11
    aj5 = -(d + y * y / rd) * d11
    if (cd ≉ zero(T))
        if (xi ≈ zero(T))
            ai4 = zero(T)
        else
            x = sqrt(xi2 + q2)
            ai4 = 1 / cdcd * (xi / rd * sdcd + 2 * atan((et * (x + q * cd) + x * (r + x) * sd) / (xi * (r + x) * cd)))
        end
        ai3 = (y * cd / rd - ale + sd * log(rd)) / cdcd
        ak1 = xi * (d11 - y11 * sd) / cd
        ak3 = (q * y11 - y * d11) / cd
        aj3 = (ak1 - aj2 * sd) / cd
        aj6 = (ak3 - aj5 * sd) / cd
    else
        rd2 = rd * rd
        ai3 = (et / rd + y * q / rd2 - ale) / 2
        ai4 = xi * y / rd2 / 2
        ak1 = xi * q / rd * d11
        ak3 = sd / rd * (xi2 * d11 - 1)
        aj3 = -xi / rd2 * (q2 * d11 - 1 / 2)
        aj6 = -y / rd2 * (xi2 * d11 - 1 / 2)
    end

    xy = xi * y11
    ai1 = -xi / rd * cd - ai4 * sd
    ai2 = log(rd) + ai3 * sd
    ak2 = 1 / r + ak3 * sd
    ak4 = xy * cd - ak1 * sd
    aj1 = aj5 * cd - aj6 * sd
    aj4 = -xy - aj2 * cd + aj3 * sd

    qx = q * x11
    qy = q * y11

    if disl[1] ≉ zero(T)
        du[1] = -xi * qy - tt - alp3 * ai1 * sd
        du[2] = -q / r + alp3 * y / rd * sd
        du[3] = q * qy - alp3 * ai2 * sd
        du[4] = xi2 * q * y32 - alp3 * aj1 * sd
        du[5] = xi * q / r3 - alp3 * aj2 * sd
        du[6] = -xi * q2 * y32 - alp3 * aj3 * sd
        du[7] = -xi * fy - d * x11 + alp3 * (xy + aj4) * sd
        du[8] = -ey + alp3 * (1 / r + aj5) * sd
        du[9] = q * fy - alp3 * (qy - aj6) * sd
        du[10] = -xi * fz - y * x11 + alp3 * ak1 * sd
        du[11] = -ez + alp3 * y * d11 * sd
        du[12] = q * fz + alp3 * ak2 * sd
        u .+= disl[1] / 2π * du
    end

    if disl[2] ≉ zero(T)
        du[1] = -q / r + alp3 * ai3 * sdcd
        du[2] = -et * qx - tt - alp3 * xi / rd * sdcd
        du[3] = q * qx + alp3 * ai4 * sdcd
        du[4] = xi * q / r3 + alp3 * aj4 * sdcd
        du[5] = et * q / r3 + qy + alp3 * aj5 * sdcd
        du[6] = -q2 / r3 + alp3 * aj6 * sdcd
        du[7] = -ey + alp3 * aj1 * sdcd
        du[8] = -et * gy - xy * sd + alp3 * aj2 * sdcd
        du[9] = q * gy + alp3 * aj3 * sdcd
        du[10] = -ez - alp3 * ak3 * sdcd
        du[11] = -et * gz - xy * cd - alp3 * xi * d11 * sdcd
        du[12] = q * gz - alp3 * ak4 * sdcd
        u .+= disl[2] / 2π * du
    end

    if disl[3] ≉ zero(T)
        du[1] = q * qy - alp3 * ai3 * sdsd
        du[2] = q * qx + alp3 * xi / rd * sdsd
        du[3] = et * qx + xi * qy - tt - alp3 * ai4 * sdsd
        du[4] = -xi * q2 * y32 - alp3 * aj4 * sdsd
        du[5] = -q2 / r3 - alp3 * aj5 * sdsd
        du[6] = q * q2 * y32 - alp3 * aj6 * sdsd
        du[7] = q * fy - alp3 * aj1 * sdsd
        du[8] = q * gy - alp3 * aj2 * sdsd
        du[9] = -q * hy - alp3 * aj3 * sdsd
        du[10] = q * fz + alp3 * ak3 * sdsd
        du[11] = q * gz + alp3 * xi * d11 * sdsd
        du[12] = -q * hz + alp3 * ak4 * sdsd
        u .+= disl[3] / 2π * du
    end
    return u
 end


function uc(sc1::B1, sc2::B2, xi::T, et::T, q::T, z::T, disl::A
    ) where {T <: Number, A <: AbstractArray{T}, B1 <: NTuple{12, T}, B2 <: NTuple{24, T}}
    alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d = sc1
    xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz = sc2

    u, du = zeros(T, 12), zeros(T, 12)

    c = d + z
    x53 = (8 * r2 + 9 * r * xi + 3 * xi2) * x11 * x11 * x11 / r2
    y53 = (8 * r2 + 9 * r * et + 3 * et2) * y11 * y11 * y11 / r2
    h = q * cd - z
    z32 = sd / r3 - h * y32
    z53 = 3 * sd / r5 - h * y53
    y0 = y11 - xi2 * y32
    z0 = z32 - xi2 * z53
    ppy = cd / r3 + q * y32 * sd
    ppz = sd / r3 - q * y32 * cd
    qq = z * y32 + z32 + z0
    qqy = 3 * c * d / r5 - qq * sd
    qqz = 3 * c * y / r5 - qq * cd + q * y32
    xy = xi * y11
    qx = q * x11
    qy = q * y11
    qr = 3 * q / r5
    cqx = c * q * x53
    cdr = (c + d) / r3
    yy0 = y / r3 - y0 * cd

    if disl[1] ≉ zero(T)
        du[1] = alp4 * xy * cd - alp5 * xi * q * z32
        du[2] = alp4 * (cd / r + 2 * qy * sd) - alp5 * c * q / r3
        du[3] = alp4 * qy * cd - alp5 * (c * et / r3 - z * y11 + xi2 * z32)
        du[4] = alp4 * y0 * cd - alp5 * q * z0
        du[5] = -alp4 * xi * (cd / r3 + 2 * q * y32 * sd) + alp5 * c * xi * qr
        du[6] = -alp4 * xi * q * y32 * cd + alp5 * xi * (3 * c * et / r5 - qq)
        du[7] = -alp4 * xi * ppy * cd - alp5 * xi * qqy
        du[8] = alp4 * 2 * (d / r3 - y0 * sd) * sd - y / r3 * cd - alp5 * (cdr * sd - et / r3 - c * y * qr)
        du[9] = -alp4 * q / r3 + yy0 * sd + alp5 * (cdr * cd + c * d * qr - (y0 * cd + q * z0) * sd)
        du[10] = alp4 * xi * ppz * cd - alp5 * xi * qqz
        du[11] = alp4 * 2 * (y / r3 - y0 * cd) * sd + d / r3 * cd - alp5 * (cdr * cd + c * d * qr)
        du[12] = yy0 * cd - alp5 * (cdr * sd - c * y * qr - y0 * sdsd + q * z0 * cd)
        u .+= disl[1] / 2π * du
    end

    if disl[2] ≉ zero(T)
        du[1] = alp4 * cd / r - qy * sd - alp5 * c * q / r3
        du[2] = alp4 * y * x11 - alp5 * c * et * q * x32
        du[3] = -d * x11 - xy * sd - alp5 * c * (x11 - q2 * x32)
        du[4] = -alp4 * xi / r3 * cd + alp5 * c * xi * qr + xi * q * y32 * sd
        du[5] = -alp4 * y / r3 + alp5 * c * et * qr
        du[6] = d / r3 - y0 * sd + alp5 * c / r3 * (1 - 3 * q2 / r2)
        du[7] = -alp4 * et / r3 + y0 * sdsd - alp5 * (cdr * sd - c * y * qr)
        du[8] = alp4 * (x11 - y * y * x32) - alp5 * c * ((d + 2 * q * cd) * x32 - y * et * q * x53)
        du[9] = xi * ppy * sd + y * d * x32 + alp5 * c * ((y + 2 * q * sd) * x32 - y * q2 * x53)
        du[10] = -q / r3 + y0 * sdcd - alp5 * (cdr * cd + c * d * qr)
        du[11] = alp4 * y * d * x32 - alp5 * c * ((y - 2 * q * sd) * x32 + d * et * q * x53)
        du[12] = -xi * ppz * sd + x11 - d * d * x32 - alp5 * c * ((d - 2 * q * cd) * x32 - d * q2 * x53)
        u .+= disl[2] / 2π * du
    end

    if disl[3] ≉ zero(T)
        du[1] = -alp4 * (sd / r + qy * cd) - alp5 * (z * y11 - q2 * z32)
        du[2] = alp4 * 2 * xy * sd + d * x11 - alp5 * c * (x11 - q2 * x32)
        du[3] = alp4 * (y * x11 + xy * cd) + alp5 * q * (c * et * x32 + xi * z32)
        du[4] = alp4 * xi / r3 * sd + xi * q * y32 * cd + alp5 * xi * (3 * c * et / r5 - 2 * z32 - z0)
        du[5] = alp4 * 2 * y0 * sd - d / r3 + alp5 * c / r3 * (1 - 3 * q2 / r2)
        du[6] = -alp4 * yy0 - alp5 * (c * et * qr - q * z0)
        du[7] = alp4 * (q / r3 + y0 * sdcd) + alp5 * (z / r3 * cd + c * d * qr - q * z0 * sd)
        du[8] = -alp4 * 2 * xi * ppy * sd - y * d * x32 + alp5 * c * ((y + 2 * q * sd) * x32 - y * q2 * x53)
        du[9] = -alp4 * (xi * ppy * cd - x11 + y * y * x32) + alp5 * (c * ((d + 2 * q * cd) * x32 - y * et * q * x53) + xi * qqy)
        du[10] = -et / r3 + y0 * cdcd - alp5 * (z / r3 * sd - c * y * qr - y0 * sdsd + q * z0 * cd)
        du[11] = alp4 * 2 * xi * ppz * sd - x11 + d * d * x32 - alp5 * c * ((d - 2 * q * cd) * x32 - d * q2 * x53)
        du[12] = alp4 * (xi * ppz * cd + y * d * x32) + alp5 * (c * ((y - 2 * q * sd) * x32 + d * et * q * x53) + xi * qqz)
        u .+= disl[3] / 2π * du
    end
    return u
end

function shared_constants_1(α::T, dip::T) where {T <: Number}
    alp1 = (1 - α) / 2
    alp2 = α / 2
    alp3 = (1 - α) / α
    alp4 = 1 - α
    alp5 = α
    sd = sind(dip)
    cd = cosd(dip)
    sdsd = sd^2
    cdcd = cd^2
    sdcd = sd * cd
    s2d = 2 * sdcd
    sc2d = cdcd - sdsd
    sc = (alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d)
end

function shared_constants_2(xi::T, et::T, q::T, sd::T, cd::T, kxi::T, ket::T) where {T <: Number}
    xi2 = xi * xi
    et2 = et * et
    q2 = q * q
    r2 = xi2 + et2 + q2
    r = sqrt(r2)

    r3 = r * r2
    r5 = r3 * r2
    y = et * cd + q * sd
    d = et * sd - q * cd

    if q ≈ zero(T)
        tt = zero(T)
    else
        tt = atan(xi * et / (q * r))
    end

    if kxi ≈ one(T)
        alx = -log(r - xi)
        x11 = zero(T)
        x32 = zero(T)
    else
        rxi = r + xi
        alx = log(rxi)
        x11 = 1 / (r*rxi)
        x32 = (r+rxi)*x11*x11/r
    end

    if ket ≈ one(T)
        ale = -log(r - et)
        y11 = zero(T)
        y32 = zero(T)
    else
        ret = r + et
        ale = log(ret)
        y11 = 1 / (r*ret)
        y32 = (r + ret) * y11 * y11/ r
    end

    ey = sd / r - y * q / r3
    ez = cd / r + d * q / r3
    fy = d / r3 + xi2 * y32 * sd
    fz = y / r3 + xi2 * y32 * cd
    gy = 2 * x11 * sd - y * q * x32
    gz = 2 * x11 * cd + d * q * x32
    hy = d * q * x32 + xi * q * y32 * sd
    hz = y * q * x32 + xi * q * y32 * cd

    return (xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz)
end
