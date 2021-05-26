# helper function for Nikkhoo's dislocation Green's func and
# Barbot's volume ones

@inline function coord_trans(x1::T, x2::T, x3::T, A::M) where {T<:Number, M<:AbstractMatrix}
    r1 = A[1,1] * x1 + A[1,2] * x2 + A[1,3] * x3
    r2 = A[2,1] * x1 + A[2,2] * x2 + A[2,3] * x3
    r3 = A[3,1] * x1 + A[3,2] * x2 + A[3,3] * x3
    r1, r2, r3
end

# @inline function coord_trans(x1::T, x2::T, x3::T, A::M) where {T<:AbstractVector, M<:AbstractMatrix}
#     r1, r2, r3 = similar(x1), similar(x2), similar(x3)
#     @inbounds @simd for i ∈ eachindex(r1)
#         r1[i] = A[1,1] * x1[i] + A[1,2] * x2[i] + A[1,3] * x3[i]
#         r2[i] = A[2,1] * x1[i] + A[2,2] * x2[i] + A[2,3] * x3[i]
#         r3[i] = A[3,1] * x1[i] + A[3,2] * x2[i] + A[3,3] * x3[i]
#     end
#     r1, r2, r3
# end

@inline function AngDisDisp(x::U, y::U, z::U, alpha::T, bx::T, by::T, bz::T, ν::T) where {T<:Number, U}
    sinA, cosA = sincos(alpha)
    eta = @. y * cosA - z * sinA
    zeta = @. y * sinA + z * cosA
    r = hypot.(x, y, z)

    if U <: Number
        zeta = ifelse(zeta > r, r, zeta)
        _z = ifelse(z > r, r, z)
    elseif U <: AbstractArray
        bool = zeta .> r
        zeta[bool] .= r[bool]
        @. bool = z > r
        _z = deepcopy(z)
        _z[bool] .= r[bool]
    end

    ux = @. bx / 8 / π / (1 - ν) * (x * y / r / (r - _z) - x * eta / r / (r - zeta))
    vx = @. bx / 8 / π / (1 - ν) * (eta * sinA / (r - zeta) - y * eta / r / (r - zeta) + y ^ 2 / r / (r - _z) + (1 - 2 * ν) * (cosA * log(r - zeta) - log(r - _z)))
    wx = @. bx / 8 / π / (1 - ν) * (eta * cosA / (r - zeta) - y / r - eta * _z / r / (r - zeta) - (1 - 2 * ν) * sinA * log(r - zeta))

    uy = @. by / 8 / π / (1 - ν) * (x ^ 2 * cosA / r / (r - zeta) - x ^ 2 / r / (r - _z) - (1 - 2 * ν) * (cosA * log(r - zeta) - log(r - _z)))
    vy = @. by * x / 8 / π / (1 - ν) * (y * cosA / r / (r - zeta) - sinA * cosA / (r - zeta) - y / r / (r - _z))
    wy = @. by * x / 8 / π / (1 - ν) * (_z * cosA / r / (r - zeta) - cosA ^ 2 / (r - zeta) + 1 / r)

    uz = @. bz * sinA / 8 / π / (1 - ν) * ((1 - 2 * ν) * log(r - zeta) - x ^ 2 / r / (r - zeta))
    vz = @. bz * x * sinA / 8 / π / (1 - ν) * (sinA / (r - zeta) - y / r / (r - zeta))
    wz = @. bz * x * sinA / 8 / π / (1 - ν) * (cosA / (r - zeta) - _z / r / (r - zeta))

    return ux + uy + uz, vx + vy + vz, wx + wy + wz
end

@inline function AngDisStrain(x::U, y::U, z::U, alpha::T, bx::T, by::T, bz::T, ν::T) where {T, U}
    sinA, cosA = sincos(alpha)
    eta = @. y * cosA - z * sinA
    zeta = @. y * sinA + z * cosA
    x2 = @. x ^ 2
    y2 = @. y ^ 2
    z2 = @. z ^ 2
    r2 = @. x2 + y2 + z2
    r = @. √(r2) # use `hypot` for higher precision
    r3 = @. r2 * r
    rz = @. r * (r - z)
    r2z2 = @. rz ^ 2
    r3z = @. r3 * (r - z)
    W = @. zeta - r
    W2 = @. W ^ 2
    Wr = @. W * r
    W2r = @. W2 * r
    Wr3 = @. W * r3
    W2r2 = @. W2 * r2
    C = @. (r * cosA - z) / Wr
    S = @. (r * sinA - y) / Wr

    rFi_rx = @. (eta / r / (r - zeta) - y / r / (r - z)) / 4 / π
    rFi_ry = @. (x / r / (r - z) - cosA * x / r / (r - zeta)) / 4 / π
    rFi_rz = @. (sinA * x / r / (r - zeta)) / 4 / π

    Exx = @. (bx * (rFi_rx) + bx / 8 / π / (1 - ν) * (eta / Wr + eta * x2 / W2r2 - eta * x2 / Wr3 + y / rz -
        x2 * y / r2z2 - x2 * y / r3z) - by * x / 8 / π / (1 - ν) * (((2 * ν + 1) / Wr + x2 / W2r2 - x2 / Wr3) * cosA +
        (2 * ν + 1) / rz - x2 / r2z2 - x2 / r3z) + bz * x * sinA / 8 / π / (1 - ν) * ((2 * ν + 1) / Wr + x2 / W2r2 - x2 / Wr3))

    Eyy = @. (by * (rFi_ry) +
        bx / 8 / π / (1 - ν) * ((1 / Wr + S ^ 2 - y2 / Wr3) * eta + (2 * ν + 1) * y / rz - y ^ 3 / r2z2 -
        y ^ 3 / r3z - 2 * ν * cosA * S) -
        by * x / 8 / π / (1 - ν) * (1 / rz - y2 / r2z2 - y2 / r3z +
        (1 / Wr + S ^ 2 - y2 / Wr3) * cosA) +
        bz * x * sinA / 8 / π / (1 - ν) * (1 / Wr + S ^ 2 - y2 / Wr3))

    Ezz = @. (bz * (rFi_rz) +
        bx / 8 / π / (1 - ν) * (eta / W / r + eta * C ^ 2 - eta * z2 / Wr3 + y * z / r3 +
        2 * ν * sinA * C) -
        by * x / 8 / π / (1 - ν) * ((1 / Wr + C ^ 2 - z2 / Wr3) * cosA + z / r3) +
        bz * x * sinA / 8 / π / (1 - ν) * (1 / Wr + C ^ 2 - z2 / Wr3))

    Exy = @. (bx * (rFi_ry) / 2 + by * (rFi_rx) / 2 -
        bx / 8 / π / (1 - ν) * (x * y2 / r2z2 - ν * x / rz + x * y2 / r3z - ν * x * cosA / Wr +
        eta * x * S / Wr + eta * x * y / Wr3) +
        by / 8 / π / (1 - ν) * (x2 * y / r2z2 - ν * y / rz + x2 * y / r3z + ν * cosA * S +
        x2 * y * cosA / Wr3 + x2 * cosA * S / Wr) -
        bz * sinA / 8 / π / (1 - ν) * (ν * S + x2 * S / Wr + x2 * y / Wr3))

    Exz = @. (bx * (rFi_rz) / 2 + bz * (rFi_rx) / 2 -
        bx / 8 / π / (1 - ν) * (-x * y / r3 + ν * x * sinA / Wr + eta * x * C / Wr +
        eta * x * z / Wr3) +
        by / 8 / π / (1 - ν) * (-x2 / r3 + ν / r + ν * cosA * C + x2 * z * cosA / Wr3 +
        x2 * cosA * C / Wr) -
        bz * sinA / 8 / π / (1 - ν) * (ν * C + x2 * C / Wr + x2 * z / Wr3))

    Eyz = @. (by * (rFi_rz) / 2 + bz * (rFi_ry) / 2 +
        bx / 8 / π / (1 - ν) * (y2 / r3 - ν / r - ν * cosA * C + ν * sinA * S + eta * sinA * cosA / W2 -
        eta * (y * cosA + z * sinA) / W2r + eta * y * z / W2r2 - eta * y * z / Wr3) -
        by * x / 8 / π / (1 - ν) * (y / r3 + sinA * cosA ^ 2 / W2 - cosA * (y * cosA + z * sinA) /
        W2r + y * z * cosA / W2r2 - y * z * cosA / Wr3) -
        bz * x * sinA / 8 / π / (1 - ν) * (y * z / Wr3 - sinA * cosA / W2 + (y * cosA + z * sinA) /
        W2r - y * z / W2r2))

    return Exx, Eyy, Ezz, Exy, Exz, Eyz
end

@inline function TensTrans(Txx1::T, Tyy1::T, Tzz1::T, Txy1::T, Txz1::T, Tyz1::T, A::M) where {T, M}
    Txx2 = @. A[1] ^ 2 * Txx1 + 2 * A[1] * A[4] * Txy1 + 2 * A[1] * A[7] * Txz1 + 2 * A[4] * A[7] * Tyz1 + A[4] ^ 2 * Tyy1 + A[7] ^ 2 * Tzz1
    Tyy2 = @. A[2] ^ 2 * Txx1 + 2 * A[2] * A[5] * Txy1 + 2 * A[2] * A[8] * Txz1 + 2 * A[5] * A[8] * Tyz1 + A[5] ^ 2 * Tyy1 + A[8] ^ 2 * Tzz1
    Tzz2 = @. A[3] ^ 2 * Txx1 + 2 * A[3] * A[6] * Txy1 + 2 * A[3] * A[9] * Txz1 + 2 * A[6] * A[9] * Tyz1 + A[6] ^ 2 * Tyy1 + A[9] ^ 2 * Tzz1
    Txy2 = @. A[1] * A[2] * Txx1 + (A[1] * A[5] + A[2] * A[4]) * Txy1 + (A[1] * A[8] + A[2] * A[7])* Txz1 + (A[8] * A[4] + A[7] * A[5]) * Tyz1 + A[5] * A[4] * Tyy1 + A[7] * A[8] * Tzz1
    Txz2 = @. A[1] * A[3] * Txx1 + (A[1] * A[6] + A[3] * A[4]) * Txy1 + (A[1] * A[9] + A[3] * A[7])* Txz1 + (A[9] * A[4] + A[7] * A[6]) * Tyz1 + A[6] * A[4] * Tyy1 + A[7] * A[9] * Tzz1
    Tyz2 = @. A[2] * A[3] * Txx1 + (A[3] * A[5] + A[2] * A[6]) * Txy1 + (A[3] * A[8] + A[2] * A[9])* Txz1 + (A[8] * A[6] + A[9] * A[5]) * Tyz1 + A[5] * A[6] * Tyy1 + A[8] * A[9] * Tzz1
    return Txx2, Tyy2, Tzz2, Txy2, Txz2, Tyz2
end

heaviside(x::T) where T = x ≤ zero(T) ? zero(T) : one(T)

xlogy(x::T, y::T) where T = isapprox(x, zero(T)) ? zero(T) : x * log(y)

xlogy(x, y) = xlogy(promote(x, y)...)

omega(x::T) where T = heaviside(x + 1/2) - heaviside(x - 1/2)

S(x::T) where T = omega(x - 1/2)
