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

dc3d_cache(T) = ntuple(x -> zeros(T, 12), Val(10))

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
function dc3d(x::T, y::T, z::T, α::T, dep::T, dip::T, al::A, aw::A, disl::AbstractVector{<:Real},
    cache=dc3d_cache(T)) where {T <: Real,A <: AbstractVecOrMat{<:Real}}

    u, du, dua, dub, duc, zerosol, xi, et, kxi, ket = cache
    fill!(u, zero(T))
    fill!(kxi, zero(T))
    fill!(ket, zero(T))

    z > zero(T) && return zerosol

    xi[1] = x - al[1]
    xi[2] = x - al[2]
    d = dep + z
    sc1 = shared_constants_1(α, dip)
    sd, cd = sc1[6], sc1[7]
    p = y * cd + d * sd
    q = y * sd - d * cd
    et[1] = p - aw[1]
    et[2] = p - aw[2]

    if q ≈ zero(T) && ((xi[1] * xi[2] ≤ zero(T) && et[1] * et[2] ≈ zero(T)) ||	(et[1] * et[2] ≤ zero(T) && xi[1] * xi[2] ≈ zero(T)))
        return zerosol
    end

    r12 = hypot(xi[1], et[2], q)
    r21 = hypot(xi[2], et[1], q)
    r22 = hypot(xi[2], et[2], q)

    (xi[1] ≤ zero(T) && (r21 + xi[2]) ≈ zero(T)) && (kxi[1] = one(T))
    (xi[1] ≤ zero(T) && (r22 + xi[2]) ≈ zero(T)) && (kxi[2] = one(T))
    (et[1] ≤ zero(T) && (r12 + et[2]) ≈ zero(T)) && (ket[1] = one(T))
    (et[1] ≤ zero(T) && (r22 + et[2]) ≈ zero(T)) && (ket[2] = one(T))

    for k in 1:2, j in 1:2
        sc2 = shared_constants_2(xi[j], et[k], q, sd, cd, kxi[k], ket[j])
        ua(dua, sc1, sc2, xi[j], et[k], q, disl)
        for i = 1:3:10
            du[i] = -dua[i]
            du[i + 1] = -dua[i + 1] * cd + dua[i + 2] * sd
            du[i + 2] = -dua[i + 1] * sd - dua[i + 2] * cd
            i < 10 && continue
            du[i] = -du[i]
            du[i + 1] = -du[i + 1]
            du[i + 2] = -du[i + 2]
        end
        for i = 1:12
            if j + k ≠ 3  u[i] += du[i] end
            if j + k == 3  u[i] -= du[i] end
        end
    end

    d = dep - z
    p = y * cd + d * sd
    q = y * sd - d * cd
    et[1] = p - aw[1]
    et[2] = p - aw[2]

    if q ≈ zero(T) && ((xi[1] * xi[2] ≤ zero(T) && et[1] * et[2] ≈ zero(T)) ||	(et[1] * et[2] ≤ zero(T) && xi[1] * xi[2] ≈ zero(T)))
        return zerosol
    end

    fill!(kxi, zero(T))
    fill!(ket, zero(T))
    r12 = hypot(xi[1], et[2], q)
    r21 = hypot(xi[2], et[1], q)
    r22 = hypot(xi[2], et[2], q)

    (xi[1] ≤ zero(T) && (r21 + xi[2]) ≈ zero(T)) && (kxi[1] = one(T))
    (xi[1] ≤ zero(T) && (r22 + xi[2]) ≈ zero(T)) && (kxi[2] = one(T))
    (et[1] ≤ zero(T) && (r12 + et[2]) ≈ zero(T)) && (ket[1] = one(T))
    (et[1] ≤ zero(T) && (r22 + et[2]) ≈ zero(T)) && (ket[2] = one(T))

    for k in 1:2, j in 1:2
        sc2 = shared_constants_2(xi[j], et[k], q, sd, cd, kxi[k], ket[j])
        ua(dua, sc1, sc2, xi[j], et[k], q, disl)
        ub(dub, sc1, sc2, xi[j], et[k], q, disl)
        uc(duc, sc1, sc2, xi[j], et[k], q, z, disl)
        for i = 1:3:10
            du[i] = dua[i] + dub[i] + z * duc[i]
            du[i + 1] = (dua[i + 1] + dub[i + 1] + z * duc[i + 1]) * cd - (dua[i + 2] + dub[i + 2] + z * duc[i + 2]) * sd
            du[i + 2] = (dua[i + 1] + dub[i + 1] - z * duc[i + 1]) * sd + (dua[i + 2] + dub[i + 2] - z * duc[i + 2]) * cd
            i < 10 && continue
            du[10] += duc[1]
            du[11] += duc[2] * cd - duc[3] * sd
            du[12] -= duc[2] * sd + duc[3] * cd
        end
        for i = 1:12
            j + k ≠ 3 && (u[i] += du[i])
            j + k == 3 && (u[i] -= du[i])
        end
    end
    return u
end

@inline function ua(u::A, sc1::B1, sc2::B2, xi::T, et::T, q::T, disl::A
    ) where {T <: Number,A <: AbstractArray{T},B1 <: NTuple{12,T},B2 <: NTuple{24,T}}
    alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d = sc1
    xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz = sc2

    fill!(u, zero(T))
    xy = xi * y11
    qx = q * x11
    qy = q * y11

    coeff = disl[1] / 2π
    u[1]  += coeff * (tt / 2 + alp2 * xi * qy)
    u[2]  += coeff * (alp2 * q / r)
    u[3]  += coeff * (alp1 * ale - alp2 * q * qy)
    u[4]  += coeff * (-alp1 * qy - alp2 * xi2 * q * y32)
    u[5]  += coeff * (-alp2 * xi * q / r3)
    u[6]  += coeff * (alp1 * xy + alp2 * xi * q2 * y32)
    u[7]  += coeff * (alp1 * xy * sd + alp2 * xi * fy + d / 2 * x11)
    u[8]  += coeff * (alp2 * ey)
    u[9]  += coeff * (alp1 * (cd / r + qy * sd) - alp2 * q * fy)
    u[10] += coeff * (alp1 * xy * cd + alp2 * xi * fz + y / 2 * x11)
    u[11] += coeff * (alp2 * ez)
    u[12] += coeff * (-alp1 * (sd / r - qy * cd) - alp2 * q * fz)

    coeff = disl[2] / 2π
    u[1]  += coeff * (alp2 * q / r)
    u[2]  += coeff * (tt / 2 + alp2 * et * qx)
    u[3]  += coeff * (alp1 * alx - alp2 * q * qx)
    u[4]  += coeff * (-alp2 * xi * q / r3)
    u[5]  += coeff * (-qy / 2 - alp2 * et * q / r3)
    u[6]  += coeff * (alp1 / r + alp2 * q2 / r3)
    u[7]  += coeff * (alp2 * ey)
    u[8]  += coeff * (alp1 * d * x11 + xy / 2 * sd + alp2 * et * gy)
    u[9]  += coeff * (alp1 * y * x11 - alp2 * q * gy)
    u[10] += coeff * (alp2 * ez)
    u[11] += coeff * (alp1 * y * x11 + xy / 2 * cd + alp2 * et * gz)
    u[12] += coeff * (-alp1 * d * x11 - alp2 * q * gz)

    coeff = disl[3] / 2π
    u[1]  += coeff * (-alp1 * ale - alp2 * q * qy)
    u[2]  += coeff * (-alp1 * alx - alp2 * q * qx)
    u[3]  += coeff * (tt / 2 - alp2 * (et * qx + xi * qy))
    u[4]  += coeff * (-alp1 * xy + alp2 * xi * q2 * y32)
    u[5]  += coeff * (-alp1 / r + alp2 * q2 / r3)
    u[6]  += coeff * (-alp1 * qy - alp2 * q * q2 * y32)
    u[7]  += coeff * (-alp1 * (cd / r + qy * sd) - alp2 * q * fy)
    u[8]  += coeff * (-alp1 * y * x11 - alp2 * q * gy)
    u[9]  += coeff * (alp1 * (d * x11 + xy * sd) + alp2 * q * hy)
    u[10] += coeff * (alp1 * (sd / r - qy * cd) - alp2 * q * fz)
    u[11] += coeff * (alp1 * d * x11 - alp2 * q * gz)
    u[12] += coeff * (alp1 * (y * x11 + xy * cd) + alp2 * q * hz)
end

@inline function ub(u::A, sc1::B1, sc2::B2, xi::T, et::T, q::T, disl::A
    ) where {T <: Number,A <: AbstractArray{T},B1 <: NTuple{12,T},B2 <: NTuple{24,T}}
    alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d = sc1
    xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz = sc2

    fill!(u, zero(T))
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

    coeff = disl[1] / 2π
    u[1]  += coeff * (-xi * qy - tt - alp3 * ai1 * sd)
    u[2]  += coeff * (-q / r + alp3 * y / rd * sd)
    u[3]  += coeff * (q * qy - alp3 * ai2 * sd)
    u[4]  += coeff * (xi2 * q * y32 - alp3 * aj1 * sd)
    u[5]  += coeff * (xi * q / r3 - alp3 * aj2 * sd)
    u[6]  += coeff * (-xi * q2 * y32 - alp3 * aj3 * sd)
    u[7]  += coeff * (-xi * fy - d * x11 + alp3 * (xy + aj4) * sd)
    u[8]  += coeff * (-ey + alp3 * (1 / r + aj5) * sd)
    u[9]  += coeff * (q * fy - alp3 * (qy - aj6) * sd)
    u[10] += coeff * (-xi * fz - y * x11 + alp3 * ak1 * sd)
    u[11] += coeff * (-ez + alp3 * y * d11 * sd)
    u[12] += coeff * (q * fz + alp3 * ak2 * sd)

    coeff = disl[2] / 2π
    u[1]  += coeff * (-q / r + alp3 * ai3 * sdcd)
    u[2]  += coeff * (-et * qx - tt - alp3 * xi / rd * sdcd)
    u[3]  += coeff * (q * qx + alp3 * ai4 * sdcd)
    u[4]  += coeff * (xi * q / r3 + alp3 * aj4 * sdcd)
    u[5]  += coeff * (et * q / r3 + qy + alp3 * aj5 * sdcd)
    u[6]  += coeff * (-q2 / r3 + alp3 * aj6 * sdcd)
    u[7]  += coeff * (-ey + alp3 * aj1 * sdcd)
    u[8]  += coeff * (-et * gy - xy * sd + alp3 * aj2 * sdcd)
    u[9]  += coeff * (q * gy + alp3 * aj3 * sdcd)
    u[10] += coeff * (-ez - alp3 * ak3 * sdcd)
    u[11] += coeff * (-et * gz - xy * cd - alp3 * xi * d11 * sdcd)
    u[12] += coeff * (q * gz - alp3 * ak4 * sdcd)

    coeff = disl[3] / 2π
    u[1]  += coeff * (q * qy - alp3 * ai3 * sdsd)
    u[2]  += coeff * (q * qx + alp3 * xi / rd * sdsd)
    u[3]  += coeff * (et * qx + xi * qy - tt - alp3 * ai4 * sdsd)
    u[4]  += coeff * (-xi * q2 * y32 - alp3 * aj4 * sdsd)
    u[5]  += coeff * (-q2 / r3 - alp3 * aj5 * sdsd)
    u[6]  += coeff * (q * q2 * y32 - alp3 * aj6 * sdsd)
    u[7]  += coeff * (q * fy - alp3 * aj1 * sdsd)
    u[8]  += coeff * (q * gy - alp3 * aj2 * sdsd)
    u[9]  += coeff * (-q * hy - alp3 * aj3 * sdsd)
    u[10] += coeff * (q * fz + alp3 * ak3 * sdsd)
    u[11] += coeff * (q * gz + alp3 * xi * d11 * sdsd)
    u[12] += coeff * (-q * hz + alp3 * ak4 * sdsd)
end


@inline function uc(u::A, sc1::B1, sc2::B2, xi::T, et::T, q::T, z::T, disl::A
    ) where {T <: Number,A <: AbstractArray{T},B1 <: NTuple{12,T},B2 <: NTuple{24,T}}
    alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d = sc1
    xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz = sc2

    fill!(u, zero(T))
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
    # qx = q * x11
    qy = q * y11
    qr = 3 * q / r5
    # cqx = c * q * x53
    cdr = (c + d) / r3
    yy0 = y / r3 - y0 * cd

    coeff = disl[1] / 2π
    u[1]  += coeff * (alp4 * xy * cd - alp5 * xi * q * z32)
    u[2]  += coeff * (alp4 * (cd / r + 2 * qy * sd) - alp5 * c * q / r3)
    u[3]  += coeff * (alp4 * qy * cd - alp5 * (c * et / r3 - z * y11 + xi2 * z32))
    u[4]  += coeff * (alp4 * y0 * cd - alp5 * q * z0)
    u[5]  += coeff * (-alp4 * xi * (cd / r3 + 2 * q * y32 * sd) + alp5 * c * xi * qr)
    u[6]  += coeff * (-alp4 * xi * q * y32 * cd + alp5 * xi * (3 * c * et / r5 - qq))
    u[7]  += coeff * (-alp4 * xi * ppy * cd - alp5 * xi * qqy)
    u[8]  += coeff * (alp4 * 2 * (d / r3 - y0 * sd) * sd - y / r3 * cd - alp5 * (cdr * sd - et / r3 - c * y * qr))
    u[9]  += coeff * (-alp4 * q / r3 + yy0 * sd + alp5 * (cdr * cd + c * d * qr - (y0 * cd + q * z0) * sd))
    u[10] += coeff * (alp4 * xi * ppz * cd - alp5 * xi * qqz)
    u[11] += coeff * (alp4 * 2 * (y / r3 - y0 * cd) * sd + d / r3 * cd - alp5 * (cdr * cd + c * d * qr))
    u[12] += coeff * (yy0 * cd - alp5 * (cdr * sd - c * y * qr - y0 * sdsd + q * z0 * cd))

    coeff = disl[2] / 2π
    u[1]  += coeff * (alp4 * cd / r - qy * sd - alp5 * c * q / r3)
    u[2]  += coeff * (alp4 * y * x11 - alp5 * c * et * q * x32)
    u[3]  += coeff * (-d * x11 - xy * sd - alp5 * c * (x11 - q2 * x32))
    u[4]  += coeff * (-alp4 * xi / r3 * cd + alp5 * c * xi * qr + xi * q * y32 * sd)
    u[5]  += coeff * (-alp4 * y / r3 + alp5 * c * et * qr)
    u[6]  += coeff * (d / r3 - y0 * sd + alp5 * c / r3 * (1 - 3 * q2 / r2))
    u[7]  += coeff * (-alp4 * et / r3 + y0 * sdsd - alp5 * (cdr * sd - c * y * qr))
    u[8]  += coeff * (alp4 * (x11 - y * y * x32) - alp5 * c * ((d + 2 * q * cd) * x32 - y * et * q * x53))
    u[9]  += coeff * (xi * ppy * sd + y * d * x32 + alp5 * c * ((y + 2 * q * sd) * x32 - y * q2 * x53))
    u[10] += coeff * (-q / r3 + y0 * sdcd - alp5 * (cdr * cd + c * d * qr))
    u[11] += coeff * (alp4 * y * d * x32 - alp5 * c * ((y - 2 * q * sd) * x32 + d * et * q * x53))
    u[12] += coeff * (-xi * ppz * sd + x11 - d * d * x32 - alp5 * c * ((d - 2 * q * cd) * x32 - d * q2 * x53))

    coeff = disl[3] / 2π
    u[1]  += coeff * (-alp4 * (sd / r + qy * cd) - alp5 * (z * y11 - q2 * z32))
    u[2]  += coeff * (alp4 * 2 * xy * sd + d * x11 - alp5 * c * (x11 - q2 * x32))
    u[3]  += coeff * (alp4 * (y * x11 + xy * cd) + alp5 * q * (c * et * x32 + xi * z32))
    u[4]  += coeff * (alp4 * xi / r3 * sd + xi * q * y32 * cd + alp5 * xi * (3 * c * et / r5 - 2 * z32 - z0))
    u[5]  += coeff * (alp4 * 2 * y0 * sd - d / r3 + alp5 * c / r3 * (1 - 3 * q2 / r2))
    u[6]  += coeff * (-alp4 * yy0 - alp5 * (c * et * qr - q * z0))
    u[7]  += coeff * (alp4 * (q / r3 + y0 * sdcd) + alp5 * (z / r3 * cd + c * d * qr - q * z0 * sd))
    u[8]  += coeff * (-alp4 * 2 * xi * ppy * sd - y * d * x32 + alp5 * c * ((y + 2 * q * sd) * x32 - y * q2 * x53))
    u[9]  += coeff * (-alp4 * (xi * ppy * cd - x11 + y * y * x32) + alp5 * (c * ((d + 2 * q * cd) * x32 - y * et * q * x53) + xi * qqy))
    u[10] += coeff * (-et / r3 + y0 * cdcd - alp5 * (z / r3 * sd - c * y * qr - y0 * sdsd + q * z0 * cd))
    u[11] += coeff * (alp4 * 2 * xi * ppz * sd - x11 + d * d * x32 - alp5 * c * ((d - 2 * q * cd) * x32 - d * q2 * x53))
    u[12] += coeff * (alp4 * (xi * ppz * cd + y * d * x32) + alp5 * (c * ((y - 2 * q * sd) * x32 + d * et * q * x53) + xi * qqz))
end

@inline function shared_constants_1(α::T, dip::T) where {T <: Number}
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
    (alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, sc2d)
end

@inline function shared_constants_2(xi::T, et::T, q::T, sd::T, cd::T, kxi::T, ket::T) where {T <: Number}
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
        x11 = 1 / (r * rxi)
        x32 = (r + rxi) * x11 * x11 / r
    end

    if ket ≈ one(T)
        ale = -log(r - et)
        y11 = zero(T)
        y32 = zero(T)
    else
        ret = r + et
        ale = log(ret)
        y11 = 1 / (r * ret)
        y32 = (r + ret) * y11 * y11 / r
    end

    ey = sd / r - y * q / r3
    ez = cd / r + d * q / r3
    fy = d / r3 + xi2 * y32 * sd
    fz = y / r3 + xi2 * y32 * cd
    gy = 2 * x11 * sd - y * q * x32
    gz = 2 * x11 * cd + d * q * x32
    hy = d * q * x32 + xi * q * y32 * sd
    hz = y * q * x32 + xi * q * y32 * cd

    (xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz)
end
