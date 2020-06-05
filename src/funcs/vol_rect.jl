#-----------------------------------------------------------------------
#  Author: James D. P. Moore (earth@jamesdpmoore.com) - 10 Jun, 2016
#  Modified: Sylvain Barbot (sbarbot@ntu.edu.sg) -
#  Earth Observatory of Singapore
#  Copyright (c) 2017 James D. P. Moore and Sylvain Barbot
#
#  This code and related code should be cited as:
#    Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
#    Associated with Distributed Anelastic Deformation in a Half Space,
#    Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
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
#-----------------------------------------------------------------------
# Translated by Pengcheng Shi (shipengcheng1230@gmail.com) 08/2019

# greens function ϵ vs u for quad4 element
export _disp_volinplane_rect, _disp_volinplane_rect!
export _strain_volinplane_rect, _strain_volinplane_rect!
export _stress_volinplane_rect, _stress_volinplane_rect!

export stress_volinplane_rect, stress_volinplane_rect!

function stress_volinplane_rect!(σ::P, x::U, z::U, qx::U, qz::U, T::U, W::U, ϕ::U, ϵxx::U, ϵxz::U, ϵzz::U, G::U, ν::U) where {U, P<:AbstractVector{<:U}}
    _stress_volinplane_rect!(σ, x, -z, qx, -qz, T, W, ϕ, ϵxx, -ϵxz, ϵzz, G, ν)
    σ[2] = -σ[2]
end

function stress_volinplane_rect(x::U, z::U, qx::U, qz::U, T::U, W::U, ϕ::U, ϵxx::U, ϵxz::U, ϵzz::U, G::U, ν::U) where {U, P}
    σ = Vector{U}(undef, 3)
    stress_volinplane_rect!(σ, x, z, qx, qz, T, W, ϕ, ϵxx, ϵxz, ϵzz, G, ν)
    return σ
end

function _strain_volinplane_rect!(ϵ::P, x2::U, x3::U, q2::U, q3::U, T::U, W::U, phi::U, epsv22p::U, epsv23p::U, epsv33p::U, G::U, nu::U) where {U, P<:AbstractVector{<:U}}
    lambda = G * 2 * nu / (1 - 2 * nu)
    epsv22 = +sind(phi) * (+epsv22p * sind(phi) + epsv23p * cosd(phi)) + cosd(phi) * (+epsv23p * sind(phi) + epsv33p * cosd(phi))
    epsv23 = +sind(phi) * (-epsv22p * cosd(phi) + epsv23p * sind(phi)) + cosd(phi) * (-epsv23p * cosd(phi) + epsv33p * sind(phi))
    epsv33 = -cosd(phi) * (-epsv22p * cosd(phi) + epsv23p * sind(phi)) + sind(phi) * (-epsv23p * cosd(phi) + epsv33p * sind(phi))
    epsvkk = epsv22 + epsv33
    x2 = x2 - q2

    let lambda=lambda, epsvkk=epsvkk, x2=x2, x3=x3, q2=q2, q3=q3, T=T, W=W, phi=phi, epsv22p=epsv22p, epsv23p=epsv23p, epsv33p=epsv33p, G=G, nu=nu

        function I223d2(y2p::R, y3p::R) where R
            (atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1) * (
              y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * (3 + (-4) * nu +
            cosd(2 * phi)) * sind(phi) + (1 / 2) * atan((-1) * (y2p + (-1) * (q3 + x3) * cosd(
              phi) + (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) *
            sind(phi)), 1) * ((-1) * (7 + 4 * nu * ((-5) + 4 * nu)) * sind(phi) + (-1) * (3 + (
            -4) * nu) * sind(3 * phi)) +(1 / 2) * ((-2) * (3 + (-4) * nu + cosd(2 * phi)) * ((
              - 1) * y2p + (q3 + (-1) * x3) * cosd(phi) + x2 * sind(phi)) * ((-1) * q3 + x3 + y2p *
            cosd(phi) + (-1) * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + (x2 + (-1) * y3p *
            cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6) * y3p + 8 *
            nu * y3p + (5 + (-8) * nu) * x2 * cosd(phi) + x2 * cosd(3 * phi) + (-2) * (q3 + (-1)
               * x3) * (4 + (-4) * nu + cosd(2 * phi)) * sind(phi) + 2 * y2p * sind(2 * phi)) +(
            q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-10) * y2p + 24 * nu *
            y2p + (-16) * nu ^ 2 * y2p + ((13 + (-28) * nu + 16 * nu ^ 2) * q3 + (11 + (-28) * nu +
            16 * nu ^ 2) * x3) * cosd(phi) + 2 * ((-3) + 4 * nu) * y2p * cosd(2 * phi) + 3 *
            q3 * cosd(3 * phi) + (-4) * nu * q3 * cosd(3 * phi) + 5 * x3 * cosd(3 * phi) + (-4)
               * nu * x3 * cosd(3 * phi) + 7 * x2 * sind(phi) + (-20) * nu * x2 * sind(phi) + 16 *
            nu ^ 2 * x2 * sind(phi) + 3 * x2 * sind(3 * phi) + (-4) * nu * x2 * sind(3 * phi)) +
            (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 *
            x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) +
            2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * (5 + 4 *
            nu * ((-3) + 2 * nu)) * y3p + (7 + 4 * nu * ((-5) + 4 * nu)) * x2 * cosd(phi) + (3 + (
            -4) * nu) * x2 * cosd(3 * phi) + (-1) * (13 * q3 + 11 * x3 + (-28) * nu * (q3 + x3) +
            16 * nu ^ 2 * (q3 + x3)) * sind(phi) + 2 * (3 + (-4) * nu) * y2p * sind(2 * phi) + ((
              - 3) * q3 + (-5) * x3 + 4 * nu * (q3 + x3)) * sind(3 * phi)) +(-2) * x3 * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
            2 * q3 * cosd(phi) + (-3) * y2p * cosd(2 * phi) + y2p * cosd(4 * phi) + (-2) * x2 *
            sind(phi) + y3p * sind(2 * phi) + 2 * x2 * sind(3 * phi) + (-1) * y3p * sind(4 *
              phi)) +(-4) * x3 * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (
            q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -2) * ((-2) * q3 * x2 * cosd(phi) + (3 * x2 * y2p + (q3 + (-1) * x3) * y3p) * cosd(
              2 * phi) + (-2) * y2p * y3p * cosd(3 * phi) + ((-1) * x2 * y2p + (q3 + x3) * y3p) *
            cosd(4 * phi) + ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (3 * q3 * y2p + x3 *
            y2p + (-1) * x2 * y3p) * sind(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2)
               * sind(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * sind(4 * phi))) +xlogy((1 / 4) * (
              (5 + (-8) * nu) * cosd(phi) + cosd(3 * phi)), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (
            -2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(
              phi) + (-2) * x3 * y3p * sind(phi)) +xlogy((1 / 4) * ((7 + 4 * nu * ((-5) + 4 * nu)
              ) * cosd(phi) + (3 + (-4) * nu) * cosd(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p * cosd(
              phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p *
            sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I223d3(y2p::R, y3p::R) where R
            ((-1) * atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1)
               * (y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * cosd(phi) * (
            3 + (-4) * nu + cosd(2 * phi)) +(1 / 2) * atan((-1) * (y2p + (-1) * (q3 + x3) *
              cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 +
            x3) * sind(phi)), 1) * ((-1) * (11 + (-28) * nu + 16 * nu ^ 2) * cosd(phi) + (-1)
               * (5 + (-4) * nu) * cosd(3 * phi)) +(1 / 4) * ((-4) * (3 + (-4) * nu + cosd(2 * phi)
              ) * ((-1) * y2p + (q3 + (-1) * x3) * cosd(phi) + x2 * sind(phi)) * ((-1) * x2 +
            y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + 2 * ((-1) * q3 +
            x3 + y2p * cosd(phi) + (-1) * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(
              phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6)
               * y3p + 8 * nu * y3p + (5 + (-8) * nu) * x2 * cosd(phi) + x2 * cosd(3 * phi) + (-2) *
            (q3 + (-1) * x3) * (4 + (-4) * nu + cosd(2 * phi)) * sind(phi) + 2 * y2p * sind(2 *
              phi)) +(-2) * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
            (-10) * y2p + 24 * nu * y2p + (-16) * nu ^ 2 * y2p + ((13 + (-28) * nu + 16 * nu ^ 2)
               * q3 + (11 + (-28) * nu + 16 * nu ^ 2) * x3) * cosd(phi) + 2 * ((-3) + 4 * nu) *
            y2p * cosd(2 * phi) + 3 * q3 * cosd(3 * phi) + (-4) * nu * q3 * cosd(3 * phi) + 5 *
            x3 * cosd(3 * phi) + (-4) * nu * x3 * cosd(3 * phi) + 7 * x2 * sind(phi) + (-20) *
            nu * x2 * sind(phi) + 16 * nu ^ 2 * x2 * sind(phi) + 3 * x2 * sind(3 * phi) + (-4) *
            nu * x2 * sind(3 * phi)) +2 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi))
               * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
            y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
               ^ (-1) * ((-2) * (5 + 4 * nu * ((-3) + 2 * nu)) * y3p + (7 + 4 * nu * ((-5) + 4 * nu)
              ) * x2 * cosd(phi) + (3 + (-4) * nu) * x2 * cosd(3 * phi) + (-1) * (13 * q3 + 11 *
            x3 + (-28) * nu * (q3 + x3) + 16 * nu ^ 2 * (q3 + x3)) * sind(phi) + 2 * (3 + (-4) *
            nu) * y2p * sind(2 * phi) + ((-3) * q3 + (-5) * x3 + 4 * nu * (q3 + x3)) * sind(3 *
              phi)) +(-8) * x3 * sind(phi) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (q3 + (-2) * y2p * cosd(phi) + 2 * (
            q3 + x3) * cosd(2 * phi) + (-1) * y2p * cosd(3 * phi) + y3p * sind(3 * phi)) +(-8)
               * x3 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * ((-2) * q3 *
            x2 * cosd(phi) + (3 * x2 * y2p + (q3 + (-1) * x3) * y3p) * cosd(2 * phi) + (-2) *
            y2p * y3p * cosd(3 * phi) + ((-1) * x2 * y2p + (q3 + x3) * y3p) * cosd(4 * phi) + ((
              - 1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (3 * q3 * y2p + x3 * y2p + (-1) * x2 *
            y3p) * sind(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2) * sind(3 * phi) +
            ((q3 + x3) * y2p + x2 * y3p) * sind(4 * phi)) +4 * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * q3 * x2 * cosd(
              phi) + (3 * x2 * y2p + (q3 + (-1) * x3) * y3p) * cosd(2 * phi) + (-2) * y2p * y3p *
            cosd(3 * phi) + ((-1) * x2 * y2p + (q3 + x3) * y3p) * cosd(4 * phi) + ((-1) *
            q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (3 * q3 * y2p + x3 * y2p + (-1) * x2 * y3p) *
            sind(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2) * sind(3 * phi) + ((q3 +
            x3) * y2p + x2 * y3p) * sind(4 * phi))) +xlogy((1 / 2) * (4 + (-4) * nu + cosd(2 *
              phi)) * sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2)
               * q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (
            -2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(
              phi)) +xlogy((1 / 4) * ((-1) * (11 + (-28) * nu + 16 * nu ^ 2) * sind(phi) + ((-5)
                +4 * nu) * sind(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (
            -2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(
              phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(
                phi)))
        end

        function I222d2(y2p::R, y3p::R) where R
            ((-1) * atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (
              y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * cosd(phi)
               * ((-3) + 4 * nu + cosd(2 * phi)) + (1 / 2) * atan((y2p + (-1) * (q3 + x3) * cosd(
                phi) + (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                  phi)) ^ (-1), 1) * ((-1) * ((-7) + 4 * (5 + (-4) * nu) * nu) * cosd(phi) + (-1) *
            (3 + (-4) * nu) * cosd(3 * phi)) +(1 / 2) * (2 * ((-3) + 4 * nu + cosd(2 * phi)) * ((
              - 1) * y3p + x2 * cosd(phi) + ((-1) * q3 + x3) * sind(phi)) * (q3 + (-1) * x3 + (-1)
               * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + (x2 + (-1) * y3p *
            cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6) * y2p + 8 *
            nu * y2p + (-2) * (q3 + (-1) * x3) * cosd(phi) * ((-4) + 4 * nu + cosd(2 * phi)) + (
            -2) * x2 * ((-2) + 4 * nu + cosd(2 * phi)) * sind(phi) + 2 * y3p * sind(2 * phi)) +(
            x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 *
            x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) +
            2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * (5 + 4 *
            nu * ((-3) + 2 * nu)) * y2p + (13 * q3 + 11 * x3 + (-28) * nu * (q3 + x3) + 16 *
            nu ^ 2 * (q3 + x3)) * cosd(phi) + ((-3) * q3 + (-5) * x3 + 4 * nu * (q3 + x3)) * cosd(
              3 * phi) + (7 + 4 * nu * ((-5) + 4 * nu)) * x2 * sind(phi) + 2 * (3 + (-4) * nu) *
            y3p * sind(2 * phi) + ((-3) + 4 * nu) * x2 * sind(3 * phi)) +(q3 + x3 + (-1) * y2p *
            cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (10 * y3p + (-24) * nu * y3p + 16 *
            nu ^ 2 * y3p + (-1) * (7 + (-20) * nu + 16 * nu ^ 2) * x2 * cosd(phi) + 2 * ((-3) +
            4 * nu) * y3p * cosd(2 * phi) + 3 * x2 * cosd(3 * phi) + (-4) * nu * x2 * cosd(3 *
              phi) + 13 * q3 * sind(phi) + (-28) * nu * q3 * sind(phi) + 16 * nu ^ 2 * q3 * sind(
                phi) + 11 * x3 * sind(phi) + (-28) * nu * x3 * sind(phi) + 16 * nu ^ 2 * x3 * sind(
                  phi) + (-3) * q3 * sind(3 * phi) + 4 * nu * q3 * sind(3 * phi) + (-5) * x3 * sind(
                    3 * phi) + 4 * nu * x3 * sind(3 * phi)) +2 * x3 * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * x2 * cosd(phi) +
            3 * y3p * cosd(2 * phi) + (-2) * x2 * cosd(3 * phi) + y3p * cosd(4 * phi) + (-2) *
            q3 * sind(phi) + y2p * sind(2 * phi) + y2p * sind(4 * phi)) +(-4) * x3 * (x2 + (-1)
               * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * ((-1) * ((-1) * q3 ^ 2 +
            x2 ^ 2 + x3 ^ 2) * cosd(phi) + ((-1) * q3 * y2p + x3 * y2p + 3 * x2 * y3p) * cosd(2 *
              phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y3p ^ 2) * cosd(3 * phi) + ((q3 + x3) *
            y2p + x2 * y3p) * cosd(4 * phi) + (-2) * q3 * x2 * sind(phi) + (x2 * y2p + (3 * q3 +
            x3) * y3p) * sind(2 * phi) + (-2) * y2p * y3p * sind(3 * phi) + (x2 * y2p + (-1) *
            (q3 + x3) * y3p) * sind(4 * phi))) +xlogy((-1 / 2) * ((-2) + 4 * nu + cosd(2 * phi)
              ) * sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) *
            q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (
            -2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(
              phi)) +xlogy((1 / 4) * ((7 + 4 * nu * ((-5) + 4 * nu)) * sind(phi) + ((-3) + 4 * nu)
               * sind(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) *
            q3 * y2p * cosd(phi) + (-2) * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) +
            (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I222d3(y2p::R, y3p::R) where R
            ((-1) * atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (
              y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * ((-3) +
            4 * nu + cosd(2 * phi)) * sind(phi) + (1 / 2) * atan((y2p + (-1) * (q3 + x3) * cosd(
              phi) + (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                phi)) ^ (-1), 1) * ((-11) * sind(phi) + (-4) * nu * ((-7) + 4 * nu) * sind(phi) +
            (-1) * ((-5) + 4 * nu) * sind(3 * phi)) +(1 / 2) * (2 * ((-3) + 4 * nu + cosd(2 *
              phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) * ((-1) *
            x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + ((-1) * q3 + x3 +
            y2p * cosd(phi) + (-1) * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(
              phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6)
               * y2p + 8 * nu * y2p + (-2) * (q3 + (-1) * x3) * cosd(phi) * ((-4) + 4 * nu + cosd(
                2 * phi)) + (-2) * x2 * ((-2) + 4 * nu + cosd(2 * phi)) * sind(phi) + 2 * y3p * sind(
                  2 * phi)) +(q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 +
            2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) *
            cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) *
            (5 + 4 * nu * ((-3) + 2 * nu)) * y2p + (13 * q3 + 11 * x3 + (-28) * nu * (q3 + x3) +
            16 * nu ^ 2 * (q3 + x3)) * cosd(phi) + ((-3) * q3 + (-5) * x3 + 4 * nu * (q3 + x3)) *
            cosd(3 * phi) + (7 + 4 * nu * ((-5) + 4 * nu)) * x2 * sind(phi) + 2 * (3 + (-4) * nu)
               * y3p * sind(2 * phi) + ((-3) + 4 * nu) * x2 * sind(3 * phi)) +(x2 + (-1) * y3p *
            cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) *
            x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-10) * y3p + 24 * nu * y3p + (
            -16) * nu ^ 2 * y3p + (7 + (-20) * nu + 16 * nu ^ 2) * x2 * cosd(phi) + (-2) * ((-3)
            +4 * nu) * y3p * cosd(2 * phi) + (-3) * x2 * cosd(3 * phi) + 4 * nu * x2 * cosd(3 *
              phi) + (-13) * q3 * sind(phi) + 28 * nu * q3 * sind(phi) + (-16) * nu ^ 2 * q3 *
            sind(phi) + (-11) * x3 * sind(phi) + 28 * nu * x3 * sind(phi) + (-16) * nu ^ 2 *
            x3 * sind(phi) + 3 * q3 * sind(3 * phi) + (-4) * nu * q3 * sind(3 * phi) + 5 * x3 *
            sind(3 * phi) + (-4) * nu * x3 * sind(3 * phi)) +4 * x3 * cosd(phi) * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
            q3 + (-2) * (q3 + x3) * cosd(2 * phi) + y2p * cosd(3 * phi) + 2 * y3p * sind(phi) + (
            -1) * y3p * sind(3 * phi)) +(-4) * x3 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p *
            sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 *
            y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) *
            sind(phi)) ^ (-2) * ((-1) * ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * cosd(phi) + ((-1)
               * q3 * y2p + x3 * y2p + 3 * x2 * y3p) * cosd(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3)
               ^ 2 + 2 * y3p ^ 2) * cosd(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * cosd(4 * phi) + (
            -2) * q3 * x2 * sind(phi) + (x2 * y2p + (3 * q3 + x3) * y3p) * sind(2 * phi) + (-2)
               * y2p * y3p * sind(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) * y3p) * sind(4 * phi))
            +2 * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
            y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
               ^ (-1) * ((-1) * ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * cosd(phi) + ((-1) * q3 * y2p +
            x3 * y2p + 3 * x2 * y3p) * cosd(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 *
            y3p ^ 2) * cosd(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * cosd(4 * phi) + (-2) *
            q3 * x2 * sind(phi) + (x2 * y2p + (3 * q3 + x3) * y3p) * sind(2 * phi) + (-2) *
            y2p * y3p * sind(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) * y3p) * sind(4 * phi))) +
            xlogy((1 / 2) * cosd(phi) * ((-4) + 4 * nu + cosd(2 * phi)), q3 ^ 2 + x2 ^ 2 + (-2) *
              q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p *
            cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 *
            y3p * sind(phi) + (-2) * x3 * y3p * sind(phi)) + xlogy((1 / 4) * ((11 + (-28) *
              nu + 16 * nu ^ 2) * cosd(phi) + ((-5) + 4 * nu) * cosd(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p *
            cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 *
            y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I233d2(y2p::R, y3p::R) where R
            (4 * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (
              -1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (
            -1), 1) * cosd(phi) + cosd(2 * phi) * (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) *
            x2 * sind(phi)) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 +
            x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 *
            y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(
              phi)) ^ (-1) + 4 * ((-1) + nu) * ((-1) + 2 * nu) * (x2 * cosd(phi) + (-1) * (q3 +
            x3) * sind(phi)) * ((-1) * q3 + (-1) * x3 + y2p * cosd(phi) + (-1) * y3p * sind(
              phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p +
            x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(
              phi)) ^ (-1) + 4 * ((-1) + nu) * ((-1) + 2 * nu) * y3p * (q3 + x3 + (-1) * y2p *
            cosd(phi) + y3p * sind(phi)) ^ (-1) * (1 + ((-1) * x2 + y3p * cosd(phi) + y2p *
            sind(phi)) ^ 2 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) ^ (-2)) ^ (
            -1) +atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1) *
              (y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * sind(phi) *
            sind(2 * phi) + ((-3) + 4 * nu) * atan((-1) * (y2p + (-1) * (q3 + x3) * cosd(phi)
              +(-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                phi)), 1) * sind(phi) * sind(2 * phi) + (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1)
               * x2 * sind(phi)) * ((-1) * q3 + x3 + y2p * cosd(phi) + (-1) * y3p * sind(phi)) *
            (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1)
               * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p)
               * sind(phi)) ^ (-1) * sind(2 * phi) + (-1) * sind(phi) * (q3 + x3 + (-1) * y2p *
            cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * ((-3) + 4 * nu) * y2p * cosd(
              phi) + ((-3) * q3 + 4 * nu * q3 + (-1) * x3 + 4 * nu * x3) * cosd(2 * phi) + ((-3) +
            4 * nu) * (q3 + (-1) * x3 + x2 * sind(2 * phi))) +(1 / 2) * (x2 + (-1) * y3p * cosd(
              phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-8) * ((-1) + nu) * ((-1) + 2 * nu)
               * y2p + (11 * q3 + x3 + 16 * nu ^ 2 * (q3 + x3) + (-4) * nu * (7 * q3 + 3 * x3)) * cosd(
              phi) + 2 * ((-3) + 4 * nu) * y2p * cosd(2 * phi) + (3 * q3 + x3 + (-4) * nu * (q3 + x3)
              ) * cosd(3 * phi) + (5 + 4 * nu * ((-5) + 4 * nu)) * x2 * sind(phi) + (3 + (-4) * nu)
               * x2 * sind(3 * phi)) +2 * x3 * sind(phi) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) *
            x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (q3 + (-2) * y2p * cosd(phi) +
            y2p * cosd(3 * phi) + 2 * x2 * sind(2 * phi) + (-1) * y3p * sind(3 * phi)) +(-2) *
            x3 * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * (((-1) *
            q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * cosd(phi) + (3 * q3 * y2p + x3 * y2p + (-1) * x2 * y3p) *
            cosd(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2) * cosd(3 * phi) + ((q3 +
            x3) * y2p + x2 * y3p) * cosd(4 * phi) + 2 * q3 * x2 * sind(phi) + ((-3) * x2 * y2p +
            ((-1) * q3 + x3) * y3p) * sind(2 * phi) + 2 * y2p * y3p * sind(3 * phi) + (x2 *
            y2p + (-1) * (q3 + x3) * y3p) * sind(4 * phi)) +xlogy((1 / 2) * cosd(2 * phi) *
              sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 *
            y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) *
            x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(phi)) +
            xlogy((1 / 4) * ((5 + 4 * nu * ((-5) + 4 * nu)) * sind(phi) + (3 + (-4) * nu) * sind(
              3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p *
            cosd(phi) + (-2) * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) *
            x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I233d3(y2p::R, y3p::R) where R
            ((-4) * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) +
              (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (
            -1), 1) * sind(phi) + atan((-1) * (y2p + (-1) * (q3 + x3) * cosd(phi) + (-1) *
              x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)), 1)
               * (3 + (-4) * nu + ((-1) + 4 * nu) * cosd(2 * phi)) * sind(phi) + cosd(2 * phi) * ((
              - 1) * y2p + (q3 + (-1) * x3) * cosd(phi) + x2 * sind(phi)) * ((-1) * q3 + x3 + y2p *
            cosd(phi) + (-1) * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + 4 * ((-1) + nu) * (
            (-1) + 2 * nu) * y3p * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) + (
            -4) * ((-1) + nu) * ((-1) + 2 * nu) * (x2 * cosd(phi) + (-1) * (q3 + x3) * sind(
              phi)) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) + (-1) * atan(
                (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1)
               * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * cosd(phi) * sind(2 * phi) +
            (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * ((-1) * x2 + y3p *
            cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 *
            y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * sind(2 * phi) + (-1) *
            sind(phi) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * ((
              - 3) + 4 * nu) * y2p * cosd(phi) + ((-3) * q3 + 4 * nu * q3 + (-1) * x3 + 4 * nu * x3)
               * cosd(2 * phi) + ((-3) + 4 * nu) * (q3 + (-1) * x3 + x2 * sind(2 * phi))) +(1 / 2) *
            (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-8) * ((-1) + nu) * ((
              - 1) + 2 * nu) * y2p + (11 * q3 + x3 + 16 * nu ^ 2 * (q3 + x3) + (-4) * nu * (7 * q3 + 3 *
            x3)) * cosd(phi) + 2 * ((-3) + 4 * nu) * y2p * cosd(2 * phi) + (3 * q3 + x3 + (-4) *
            nu * (q3 + x3)) * cosd(3 * phi) + (5 + 4 * nu * ((-5) + 4 * nu)) * x2 * sind(phi) + (
            3 + (-4) * nu) * x2 * sind(3 * phi)) +2 * x3 * cosd(phi) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (q3 + 2 * x3 + (
            -2) * (q3 + x3) * cosd(2 * phi) + y2p * cosd(3 * phi) + 2 * y3p * sind(phi) + (-1) *
            y3p * sind(3 * phi)) +(-2) * x3 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(
              phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p +
            x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(
              phi)) ^ (-2) * (((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * cosd(phi) + (3 * q3 * y2p + x3 *
            y2p + (-1) * x2 * y3p) * cosd(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2)
               * cosd(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * cosd(4 * phi) + 2 * q3 * x2 * sind(
              phi) + ((-3) * x2 * y2p + ((-1) * q3 + x3) * y3p) * sind(2 * phi) + 2 * y2p * y3p *
            sind(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) * y3p) * sind(4 * phi)) +(q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
            ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * cosd(phi) + (3 * q3 * y2p + x3 * y2p + (-1) * x2 *
            y3p) * cosd(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2) * cosd(3 * phi) +
            ((q3 + x3) * y2p + x2 * y3p) * cosd(4 * phi) + 2 * q3 * x2 * sind(phi) + ((-3) *
            x2 * y2p + ((-1) * q3 + x3) * y3p) * sind(2 * phi) + 2 * y2p * y3p * sind(3 * phi) +
            (x2 * y2p + (-1) * (q3 + x3) * y3p) * sind(4 * phi)) +xlogy((-1 / 2) * cosd(phi)
               * cosd(2 * phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) *
            q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (
            -2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(
              phi)) +xlogy((1 / 4) * ((1 + (-12) * nu + 16 * nu ^ 2) * cosd(phi) + (1 + (-4) * nu)
               * cosd(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) *
            q3 * y2p * cosd(phi) + (-2) * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) +
            (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I232d2(y2p::R, y3p::R) where R
            ((-4) * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) +
              (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                phi)), 1) * sind(phi) + cosd(2 * phi) * (y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) *
            x3) * sind(phi)) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 +
            x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 *
            y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(
              phi)) ^ (-1) + (-4) * ((-1) + nu) * ((-1) + 2 * nu) * ((q3 + x3) * cosd(phi) + x2 *
            sind(phi)) * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) +
            4 * ((-1) + nu) * ((-1) + 2 * nu) * y2p * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p *
            sind(phi)) ^ (-1) * (1 + ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) ^ 2 * (
            q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) ^ (-2)) ^ (-1) + (-1) * atan(
              (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 *
            cosd(phi) + (q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * cosd(phi) * sind(2 * phi) +
            ((-3) + 4 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (-1) * x2 * sind(
              phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (-1), 1) * cosd(
                phi) * sind(2 * phi) + ((-1) * y3p + x2 * cosd(phi) + ((-1) * q3 + x3) * sind(phi))
               * (q3 + (-1) * x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (
            -2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 *
            y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (
            -1) * sind(2 * phi) + cosd(phi) * ((-1) * q3 + (-1) * x3 + y2p * cosd(phi) + (-1) *
            y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (
            q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p)
               * sind(phi)) ^ (-1) * ((((-3) + 4 * nu) * q3 + ((-1) + 4 * nu) * x3) * cosd(2 *
              phi) + (-1) * ((-3) + 4 * nu) * (q3 + (-1) * x3 + 2 * y3p * sind(phi) + (-1) * x2 *
            sind(2 * phi))) +(1 / 2) * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi))
               * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
            y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
               ^ (-1) * (8 * ((-1) + nu) * ((-1) + 2 * nu) * y3p + ((-5) + 4 * (5 + (-4) * nu) *
            nu) * x2 * cosd(phi) + 2 * ((-3) + 4 * nu) * y3p * cosd(2 * phi) + (3 + (-4) * nu) *
            x2 * cosd(3 * phi) + (11 * q3 + x3 + 16 * nu ^ 2 * (q3 + x3) + (-4) * nu * (7 * q3 + 3 *
            x3)) * sind(phi) + ((-3) * q3 + (-1) * x3 + 4 * nu * (q3 + x3)) * sind(3 * phi)) +(
            -2) * x3 * cosd(phi) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2)
               * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) *
            y3p) * sind(phi)) ^ (-1) * (q3 + (-1) * y2p * cosd(3 * phi) + 2 * y3p * sind(phi)
            +(-2) * x2 * sind(2 * phi) + y3p * sind(3 * phi)) +(-2) * x3 * (x2 + (-1) * y3p *
            cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) *
            x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * ((-2) * q3 * x2 * cosd(phi) + (
            x2 * y2p + (3 * q3 + x3) * y3p) * cosd(2 * phi) + (-2) * y2p * y3p * cosd(3 * phi) +
            (x2 * y2p + (-1) * (q3 + x3) * y3p) * cosd(4 * phi) + ((-1) * q3 ^ 2 + x2 ^ 2 +
            x3 ^ 2) * sind(phi) + (q3 * y2p + (-1) * x3 * y2p + (-3) * x2 * y3p) * sind(2 *
              phi) + (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y3p ^ 2) * sind(3 * phi) + (-1) * ((q3 + x3) *
            y2p + x2 * y3p) * sind(4 * phi)) +xlogy((1 / 2) * cosd(phi) * cosd(2 * phi),
              q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(
                phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p *
            sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(phi)) + xlogy((1 / 4)
               * (((-5) + 4 * (5 + (-4) * nu) * nu) * cosd(phi) + (3 + (-4) * nu) * cosd(3 * phi))
              , q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi)
            +(-2) * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p *
            sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I232d3(y2p::R, y3p::R) where R
            ((-4) * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) +
              (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                phi)), 1) * cosd(phi) + (-1) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (-1) *
                  x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (-1), 1)
               * cosd(phi) * (3 + (-4) * nu + (1 + (-4) * nu) * cosd(2 * phi)) + cosd(2 * phi) * (
            y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) * (q3 + (-1) * x3 + (-1)
               * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + 4 * ((-1) + nu) * (
            (-1) + 2 * nu) * y2p * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) + (
            -4) * ((-1) + nu) * ((-1) + 2 * nu) * ((q3 + x3) * cosd(phi) + x2 * sind(phi)) * ((
              - 1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) + (-1) * atan((y2p + ((
              - 1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi)
            +(q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * sind(phi) * sind(2 * phi) + (y3p + (
            -1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) * ((-1) * x2 + y3p * cosd(
              phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 *
            y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * sind(2 * phi) + (-1) *
            cosd(phi) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((((-3) + 4 *
            nu) * q3 + ((-1) + 4 * nu) * x3) * cosd(2 * phi) + (-1) * ((-3) + 4 * nu) * (q3 + (
            -1) * x3 + 2 * y3p * sind(phi) + (-1) * x2 * sind(2 * phi))) +(1 / 2) * (q3 + x3 + (
            -1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) *
            x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (8 * ((-1) + nu) * ((-1) + 2 *
            nu) * y3p + ((-5) + 4 * (5 + (-4) * nu) * nu) * x2 * cosd(phi) + 2 * ((-3) + 4 * nu)
               * y3p * cosd(2 * phi) + (3 + (-4) * nu) * x2 * cosd(3 * phi) + (11 * q3 + x3 + 16 *
            nu ^ 2 * (q3 + x3) + (-4) * nu * (7 * q3 + 3 * x3)) * sind(phi) + ((-3) * q3 + (-1) *
            x3 + 4 * nu * (q3 + x3)) * sind(3 * phi)) +2 * x3 * sind(phi) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (q3 + 2 * x3 + (
            -2) * y2p * cosd(phi) + 2 * (q3 + x3) * cosd(2 * phi) + (-1) * y2p * cosd(3 * phi) +
            y3p * sind(3 * phi)) +(-2) * x3 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(
              phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p +
            x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(
              phi)) ^ (-2) * ((-2) * q3 * x2 * cosd(phi) + (x2 * y2p + (3 * q3 + x3) * y3p) *
            cosd(2 * phi) + (-2) * y2p * y3p * cosd(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) *
            y3p) * cosd(4 * phi) + ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (q3 * y2p + (
            -1) * x3 * y2p + (-3) * x2 * y3p) * sind(2 * phi) + (x2 ^ 2 + (q3 + x3) ^ 2 + 2 *
            y3p ^ 2) * sind(3 * phi) + (-1) * ((q3 + x3) * y2p + x2 * y3p) * sind(4 * phi)) +(
            q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -1) * ((-2) * q3 * x2 * cosd(phi) + (x2 * y2p + (3 * q3 + x3) * y3p) * cosd(2 *
              phi) + (-2) * y2p * y3p * cosd(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) * y3p) *
            cosd(4 * phi) + ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (q3 * y2p + (-1) *
            x3 * y2p + (-3) * x2 * y3p) * sind(2 * phi) + (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y3p ^ 2) *
            sind(3 * phi) + (-1) * ((q3 + x3) * y2p + x2 * y3p) * sind(4 * phi)) +xlogy((1 / 2)
               * cosd(2 * phi) * sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 *
            y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) *
            x3 * y3p * sind(phi)) + xlogy((1 / 4) * ((1 + (-12) * nu + 16 * nu ^ 2) * sind(phi)
              +((-1) + 4 * nu) * sind(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p * cosd(phi) + (-2) * x2 *
            y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 *
            y3p * sind(phi)))
        end

        function I323d2(y2p::R, y3p::R) where R
            ((-4) * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) +
              (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (
            -1), 1) * cosd(phi) + cosd(2 * phi) * (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) *
            x2 * sind(phi)) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 +
            x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 *
            y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(
              phi)) ^ (-1) + (-4) * ((-1) + nu) * ((-1) + 2 * nu) * (x2 * cosd(phi) + (-1) * (
            q3 + x3) * sind(phi)) * ((-1) * q3 + (-1) * x3 + y2p * cosd(phi) + (-1) * y3p *
            sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 *
            y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) *
            sind(phi)) ^ (-1) + (-4) * ((-1) + nu) * ((-1) + 2 * nu) * y3p * (q3 + x3 + (-1) *
            y2p * cosd(phi) + y3p * sind(phi)) ^ (-1) * (1 + ((-1) * x2 + y3p * cosd(phi) +
            y2p * sind(phi)) ^ 2 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) ^ (-2)
              ) ^ (-1) + atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (
              -1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * sind(phi)
               * sind(2 * phi) + ((-3) + 4 * nu) * atan((-1) * (y2p + (-1) * (q3 + x3) * cosd(
                phi) + (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) *
              sind(phi)), 1) * sind(phi) * sind(2 * phi) + (y2p + ((-1) * q3 + x3) * cosd(phi) + (
            -1) * x2 * sind(phi)) * ((-1) * q3 + x3 + y2p * cosd(phi) + (-1) * y3p * sind(phi)
              ) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (
            -1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 *
            y3p) * sind(phi)) ^ (-1) * sind(2 * phi) + (-1) * sind(phi) * (q3 + x3 + (-1) *
            y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * ((-3) + 4 * nu) * y2p * cosd(
              phi) + ((-3) * q3 + 4 * nu * q3 + (-5) * x3 + 4 * nu * x3) * cosd(2 * phi) + ((-3) +
            4 * nu) * (q3 + (-1) * x3 + x2 * sind(2 * phi))) +(1 / 2) * (x2 + (-1) * y3p * cosd(
              phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (8 * ((-1) + nu) * ((-1) + 2 * nu) *
            y2p + (-1) * ((5 + 4 * nu * ((-5) + 4 * nu)) * q3 + (19 + 4 * nu * ((-9) + 4 * nu)) *
            x3) * cosd(phi) + 2 * ((-3) + 4 * nu) * y2p * cosd(2 * phi) + (3 * q3 + 5 * x3 + (-4)
               * nu * (q3 + x3)) * cosd(3 * phi) + ((-11) + 4 * (7 + (-4) * nu) * nu) * x2 * sind(
              phi) + (3 + (-4) * nu) * x2 * sind(3 * phi)) +(-2) * x3 * sind(phi) * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
            q3 + (-2) * y2p * cosd(phi) + y2p * cosd(3 * phi) + 2 * x2 * sind(2 * phi) + (-1) *
            y3p * sind(3 * phi)) +(-2) * x3 * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p *
            sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 *
            y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) *
            sind(phi)) ^ (-2) * ((q3 ^ 2 + (-1) * x2 ^ 2 + (-1) * x3 ^ 2) * cosd(phi) + ((-1)
               * (3 * q3 + x3) * y2p + x2 * y3p) * cosd(2 * phi) + (x2 ^ 2 + (q3 + x3) ^ 2 + 2 *
            y2p ^ 2) * cosd(3 * phi) + (-1) * ((q3 + x3) * y2p + x2 * y3p) * cosd(4 * phi) + (
            -2) * q3 * x2 * sind(phi) + (3 * x2 * y2p + (q3 + (-1) * x3) * y3p) * sind(2 * phi)
            +(-2) * y2p * y3p * sind(3 * phi) + ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(4 *
              phi)) +xlogy((1 / 2) * cosd(2 * phi) * sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 +
                x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (
              -2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(
                phi) + (-2) * x3 * y3p * sind(phi)) +xlogy((1 / 4) * (((-11) + 4 * (7 + (-4) * nu)
               * nu) * sind(phi) + (3 + (-4) * nu) * sind(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p * cosd(
              phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p *
            sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I323d3(y2p::R, y3p::R) where R
            (4 * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (
              -1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (
            -1), 1) * sind(phi) + atan((-1) * (y2p + (-1) * (q3 + x3) * cosd(phi) + (-1) *
              x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)), 1)
               * (3 + (-4) * nu + ((-5) + 4 * nu) * cosd(2 * phi)) * sind(phi) + cosd(2 * phi) * ((
              - 1) * y2p + (q3 + (-1) * x3) * cosd(phi) + x2 * sind(phi)) * ((-1) * q3 + x3 + y2p *
            cosd(phi) + (-1) * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + (-4) * ((-1) + nu)
               * ((-1) + 2 * nu) * y3p * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (
            q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -1) +4 * ((-1) + nu) * ((-1) + 2 * nu) * (x2 * cosd(phi) + (-1) * (q3 + x3) * sind(
              phi)) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) + (-1) * atan(
                (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1)
               * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * cosd(phi) * sind(2 * phi) +
            (y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * ((-1) * x2 + y3p *
            cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 *
            y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * sind(2 * phi) + (-1) *
            sind(phi) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * ((
              - 3) + 4 * nu) * y2p * cosd(phi) + ((-3) * q3 + 4 * nu * q3 + (-5) * x3 + 4 * nu * x3)
               * cosd(2 * phi) + ((-3) + 4 * nu) * (q3 + (-1) * x3 + x2 * sind(2 * phi))) +(1 / 2) *
            (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (8 * ((-1) + nu) * ((-1)
            +2 * nu) * y2p + (-1) * ((5 + 4 * nu * ((-5) + 4 * nu)) * q3 + (19 + 4 * nu * ((-9) +
            4 * nu)) * x3) * cosd(phi) + 2 * ((-3) + 4 * nu) * y2p * cosd(2 * phi) + (3 * q3 +
            5 * x3 + (-4) * nu * (q3 + x3)) * cosd(3 * phi) + ((-11) + 4 * (7 + (-4) * nu) * nu)
               * x2 * sind(phi) + (3 + (-4) * nu) * x2 * sind(3 * phi)) +(-2) * x3 * cosd(phi) *
            (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -1) * (q3 + 2 * x3 + (-2) * (q3 + x3) * cosd(2 * phi) + y2p * cosd(3 * phi) + 2 *
            y3p * sind(phi) + (-1) * y3p * sind(3 * phi)) +(-2) * x3 * (q3 + x3 + (-1) * y2p *
            cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * ((q3 ^ 2 + (-1) * x2 ^ 2 + (-1) *
            x3 ^ 2) * cosd(phi) + ((-1) * (3 * q3 + x3) * y2p + x2 * y3p) * cosd(2 * phi) + (
            x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2) * cosd(3 * phi) + (-1) * ((q3 + x3) * y2p + x2 *
            y3p) * cosd(4 * phi) + (-2) * q3 * x2 * sind(phi) + (3 * x2 * y2p + (q3 + (-1) * x3)
               * y3p) * sind(2 * phi) + (-2) * y2p * y3p * sind(3 * phi) + ((-1) * x2 * y2p + (
            q3 + x3) * y3p) * sind(4 * phi)) +(q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((q3 ^ 2 + (-1) * x2 ^ 2 + (-1) *
            x3 ^ 2) * cosd(phi) + ((-1) * (3 * q3 + x3) * y2p + x2 * y3p) * cosd(2 * phi) + (
            x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2) * cosd(3 * phi) + (-1) * ((q3 + x3) * y2p + x2 *
            y3p) * cosd(4 * phi) + (-2) * q3 * x2 * sind(phi) + (3 * x2 * y2p + (q3 + (-1) * x3)
               * y3p) * sind(2 * phi) + (-2) * y2p * y3p * sind(3 * phi) + ((-1) * x2 * y2p + (
            q3 + x3) * y3p) * sind(4 * phi)) +xlogy((-1 / 2) * cosd(phi) * cosd(2 * phi),
              q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(
                phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p *
            sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(phi)) + xlogy((1 / 4)
               * ((-1) * (19 + 4 * nu * ((-9) + 4 * nu)) * cosd(phi) + (5 + (-4) * nu) * cosd(3 *
              phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p *
            cosd(phi) + (-2) * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) *
            x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I322d2(y2p::R, y3p::R) where R
            (4 * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (
              -1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                phi)), 1) * sind(phi) + 2 * ((-3) + 4 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(
                  phi) + (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                    phi)) ^ (-1), 1) * cosd(phi) ^ 2 * sind(phi) + cosd(2 * phi) * (y3p + (-1) * x2 *
            cosd(phi) + (q3 + (-1) * x3) * sind(phi)) * ((-1) * x2 + y3p * cosd(phi) + y2p *
            sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (
            q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 *
            y3p + x3 * y3p) * sind(phi)) ^ (-1) + 4 * ((-1) + nu) * ((-1) + 2 * nu) * ((q3 + x3)
               * cosd(phi) + x2 * sind(phi)) * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)
              ) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
            y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
               ^ (-1) + (-4) * ((-1) + nu) * ((-1) + 2 * nu) * y2p * (q3 + x3 + (-1) * y2p * cosd(
              phi) + y3p * sind(phi)) ^ (-1) * (1 + ((-1) * x2 + y3p * cosd(phi) + y2p * sind(
                phi)) ^ 2 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) ^ (-2)) ^ (-1) + (
            -1) * atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (y3p +
              (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * cosd(phi) *
            sind(2 * phi) + ((-1) * y3p + x2 * cosd(phi) + ((-1) * q3 + x3) * sind(phi)) * (q3 +
            (-1) * x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) *
            cosd(phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) *
            sind(2 * phi) + (-1) * cosd(phi) * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(
              phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p +
            x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(
              phi)) ^ (-1) * ((((-3) + 4 * nu) * q3 + ((-5) + 4 * nu) * x3) * cosd(2 * phi) + (
            -1) * ((-3) + 4 * nu) * (q3 + (-1) * x3 + 2 * y3p * sind(phi) + (-1) * x2 * sind(2 *
              phi))) +(1 / 2) * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (
            q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -1) * ((-8) * ((-1) + nu) * ((-1) + 2 * nu) * y3p + (11 + 4 * nu * ((-7) + 4 * nu))
               * x2 * cosd(phi) + 2 * ((-3) + 4 * nu) * y3p * cosd(2 * phi) + (3 + (-4) * nu) *
            x2 * cosd(3 * phi) + (-1) * (5 * q3 + 19 * x3 + 16 * nu ^ 2 * (q3 + x3) + (-4) * nu * (
            5 * q3 + 9 * x3)) * sind(phi) + ((-3) * q3 + (-5) * x3 + 4 * nu * (q3 + x3)) * sind(
              3 * phi)) +2 * x3 * cosd(phi) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (q3 + (-1) * y2p * cosd(3 * phi) + 2 *
            y3p * sind(phi) + (-2) * x2 * sind(2 * phi) + y3p * sind(3 * phi)) +(-2) * x3 * (
            x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 *
            x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) +
            2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * (2 * q3 * x2 * cosd(
              phi) + (-1) * (x2 * y2p + (3 * q3 + x3) * y3p) * cosd(2 * phi) + 2 * y2p * y3p *
            cosd(3 * phi) + ((-1) * x2 * y2p + (q3 + x3) * y3p) * cosd(4 * phi) + (-1) * ((-1)
               * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + ((-1) * q3 * y2p + x3 * y2p + 3 * x2 * y3p)
               * sind(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y3p ^ 2) * sind(3 * phi) + ((
              q3+ x3) * y2p + x2 * y3p) * sind(4 * phi)) +xlogy((1 / 2) * cosd(phi) * cosd(2 *
                phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p *
              cosd(phi) + 2 * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 *
              y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(phi)) + xlogy(
                (1 / 4) * ((11 + 4 * nu * ((-7) + 4 * nu)) * cosd(phi) + (3 + (-4) * nu) * cosd(3 *
                  phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p *
              cosd(phi) + (-2) * x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) *
              x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I322d3(y2p::R, y3p::R) where R
            (4 * ((-1) + nu) * ((-1) + 2 * nu) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (
              -1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(
                phi)), 1) * cosd(phi) + atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (-1) * x2 *
                  sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (-1), 1) *
            cosd(phi) * ((-3) + 4 * nu + ((-5) + 4 * nu) * cosd(2 * phi)) + cosd(2 * phi) * (
            y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) * (q3 + (-1) * x3 + (-1)
               * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
               * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + (-4) * ((-1) + nu)
               * ((-1) + 2 * nu) * y2p * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (
            q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -1) +4 * ((-1) + nu) * ((-1) + 2 * nu) * ((q3 + x3) * cosd(phi) + x2 * sind(phi)) *
            ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) + (-1) * atan((y2p + ((
              - 1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi)
            +(q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * sind(phi) * sind(2 * phi) + (y3p + (
            -1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) * ((-1) * x2 + y3p * cosd(
              phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 *
            y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * sind(2 * phi) + (-1) *
            cosd(phi) * ((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(
              phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((((-3) + 4 *
            nu) * q3 + ((-5) + 4 * nu) * x3) * cosd(2 * phi) + (-1) * ((-3) + 4 * nu) * (q3 + (
            -1) * x3 + 2 * y3p * sind(phi) + (-1) * x2 * sind(2 * phi))) +(1 / 2) * (q3 + x3 + (
            -1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) *
            x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-8) * ((-1) + nu) * ((-1) +
            2 * nu) * y3p + (11 + 4 * nu * ((-7) + 4 * nu)) * x2 * cosd(phi) + 2 * ((-3) + 4 * nu)
               * y3p * cosd(2 * phi) + (3 + (-4) * nu) * x2 * cosd(3 * phi) + (-1) * (5 * q3 + 19 *
            x3 + 16 * nu ^ 2 * (q3 + x3) + (-4) * nu * (5 * q3 + 9 * x3)) * sind(phi) + ((-3) *
            q3 + (-5) * x3 + 4 * nu * (q3 + x3)) * sind(3 * phi)) +(-2) * x3 * sind(phi) * (
            q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -1) * (q3 + 2 * x3 + (-2) * y2p * cosd(phi) + 2 * (q3 + x3) * cosd(2 * phi) + (-1) *
            y2p * cosd(3 * phi) + y3p * sind(3 * phi)) +(-2) * x3 * (q3 + x3 + (-1) * y2p *
            cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * (2 * q3 * x2 * cosd(phi) + (-1) * (
            x2 * y2p + (3 * q3 + x3) * y3p) * cosd(2 * phi) + 2 * y2p * y3p * cosd(3 * phi) + ((
              - 1) * x2 * y2p + (q3 + x3) * y3p) * cosd(4 * phi) + (-1) * ((-1) * q3 ^ 2 + x2 ^ 2 +
            x3 ^ 2) * sind(phi) + ((-1) * q3 * y2p + x3 * y2p + 3 * x2 * y3p) * sind(2 * phi) + (
            -1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y3p ^ 2) * sind(3 * phi) + ((q3 + x3) * y2p + x2 *
            y3p) * sind(4 * phi)) +(q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2)
               * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) *
            y3p) * sind(phi)) ^ (-1) * (2 * q3 * x2 * cosd(phi) + (-1) * (x2 * y2p + (3 * q3 +
            x3) * y3p) * cosd(2 * phi) + 2 * y2p * y3p * cosd(3 * phi) + ((-1) * x2 * y2p + (
            q3 + x3) * y3p) * cosd(4 * phi) + (-1) * ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(
              phi) + ((-1) * q3 * y2p + x3 * y2p + 3 * x2 * y3p) * sind(2 * phi) + (-1) * (x2 ^ 2 +
            (q3 + x3) ^ 2 + 2 * y3p ^ 2) * sind(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * sind(4 *
              phi)) +xlogy((1 / 2) * cosd(2 * phi) * sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 +
                x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(phi) + (
              -2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p * sind(
                phi) + (-2) * x3 * y3p * sind(phi)) +xlogy((1 / 4) * ((-1) * (19 + (-36) * nu +
                  16 * nu ^ 2) * sind(phi) + ((-5) + 4 * nu) * sind(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p *
            cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 *
            y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I333d2(y2p::R, y3p::R) where R
            ((-1) * atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1)
               * (y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * ((-3) + 4 *
            nu + cosd(2 * phi)) * sind(phi) + (1 / 2) * atan((-1) * (y2p + (-1) * (q3 + x3) *
              cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 +
            x3) * sind(phi)), 1) * (((-13) + 4 * (7 + (-4) * nu) * nu) * sind(phi) + (3 + (-4)
               * nu) * sind(3 * phi)) +(1 / 2) * (2 * ((-3) + 4 * nu + cosd(2 * phi)) * ((-1) *
            y2p + (q3 + (-1) * x3) * cosd(phi) + x2 * sind(phi)) * ((-1) * q3 + x3 + y2p * cosd(
              phi) + (-1) * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 *
            y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + (x2 + (-1) * y3p * cosd(
              phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 *
            y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6) * y3p + 8 * nu *
            y3p + (7 + (-8) * nu) * x2 * cosd(phi) + (-1) * x2 * cosd(3 * phi) + 2 * (q3 + (-1) *
            x3) * ((-2) + 4 * nu + cosd(2 * phi)) * sind(phi) + (-2) * y2p * sind(2 * phi)) +(
            q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-10) * y2p + 24 * nu *
            y2p + (-16) * nu ^ 2 * y2p + ((7 + (-20) * nu + 16 * nu ^ 2) * q3 + (5 + (-20) * nu +
            16 * nu ^ 2) * x3) * cosd(phi) + (-2) * ((-3) + 4 * nu) * y2p * cosd(2 * phi) + (
            -3) * q3 * cosd(3 * phi) + 4 * nu * q3 * cosd(3 * phi) + (-1) * x3 * cosd(3 * phi) +
            4 * nu * x3 * cosd(3 * phi) + 13 * x2 * sind(phi) + (-28) * nu * x2 * sind(phi) +
            16 * nu ^ 2 * x2 * sind(phi) + (-3) * x2 * sind(3 * phi) + 4 * nu * x2 * sind(3 *
              phi)) +(x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 +
            2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) *
            cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) *
            (5 + 4 * nu * ((-3) + 2 * nu)) * y3p + (13 + 4 * nu * ((-7) + 4 * nu)) * x2 * cosd(
              phi) + ((-3) + 4 * nu) * x2 * cosd(3 * phi) + (-1) * (7 * q3 + 5 * x3 + (-20) * nu * (
            q3 + x3) +16 * nu ^ 2 * (q3 + x3)) * sind(phi) + 2 * ((-3) + 4 * nu) * y2p * sind(2 *
              phi) + (3 * q3 + x3 + (-4) * nu * (q3 + x3)) * sind(3 * phi)) +(-2) * x3 * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
            2 * q3 * cosd(phi) + (-3) * y2p * cosd(2 * phi) + y2p * cosd(4 * phi) + (-2) * x2 *
            sind(phi) + y3p * sind(2 * phi) + 2 * x2 * sind(3 * phi) + (-1) * y3p * sind(4 *
              phi)) +(-4) * x3 * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (
            q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p +
            x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (
            -2) * ((-2) * q3 * x2 * cosd(phi) + (3 * x2 * y2p + (q3 + (-1) * x3) * y3p) * cosd(
              2 * phi) + (-2) * y2p * y3p * cosd(3 * phi) + ((-1) * x2 * y2p + (q3 + x3) * y3p) *
            cosd(4 * phi) + ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (3 * q3 * y2p + x3 *
            y2p + (-1) * x2 * y3p) * sind(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y2p ^ 2)
               * sind(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * sind(4 * phi))) +xlogy((1 / 4) * (
              (7 + (-8) * nu) * cosd(phi) + (-1) * cosd(3 * phi)), q3 ^ 2 + x2 ^ 2 + (-2) * q3 *
            x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(
              phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p *
            sind(phi) + (-2) * x3 * y3p * sind(phi)) +xlogy((1 / 4) * ((13 + 4 * nu * ((-7) +
              4 * nu)) * cosd(phi) + ((-3) + 4 * nu) * cosd(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 * q3 *
            x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p * cosd(
              phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p *
            sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I333d3(y2p::R, y3p::R) where R
            (atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) ^ (-1) * (
              y3p + (-1) * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)), 1) * cosd(phi) * ((
              - 3) +4 * nu + cosd(2 * phi)) +(1 / 2) * atan((-1) * (y2p + (-1) * (q3 + x3) * cosd(
                phi) + (-1) * x2 * sind(phi)) ^ (-1) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) *
              sind(phi)), 1) * ((-1) * (5 + (-20) * nu + 16 * nu ^ 2) * cosd(phi) + (1 + (-4) *
            nu) * cosd(3 * phi)) +(1 / 2) * (2 * ((-3) + 4 * nu + cosd(2 * phi)) * ((-1) * y2p +
            (q3 + (-1) * x3) * cosd(phi) + x2 * sind(phi)) * ((-1) * x2 + y3p * cosd(phi) +
            y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2)
               * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1) *
            q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + ((-1) * q3 + x3 + y2p * cosd(phi) + (-1)
               * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (
            -2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1)
               * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6) * y3p + 8 * nu * y3p + (7 + (-8)
               * nu) * x2 * cosd(phi) + (-1) * x2 * cosd(3 * phi) + 2 * (q3 + (-1) * x3) * ((-2) +
            4 * nu + cosd(2 * phi)) * sind(phi) + (-2) * y2p * sind(2 * phi)) +((-1) * x2 +
            y3p * cosd(phi) + y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-10) * y2p + 24 * nu * y2p + (-16)
               * nu ^ 2 * y2p + ((7 + (-20) * nu + 16 * nu ^ 2) * q3 + (5 + (-20) * nu + 16 * nu ^ 2)
               * x3) * cosd(phi) + (-2) * ((-3) + 4 * nu) * y2p * cosd(2 * phi) + (-3) * q3 *
            cosd(3 * phi) + 4 * nu * q3 * cosd(3 * phi) + (-1) * x3 * cosd(3 * phi) + 4 * nu *
            x3 * cosd(3 * phi) + 13 * x2 * sind(phi) + (-28) * nu * x2 * sind(phi) + 16 *
            nu ^ 2 * x2 * sind(phi) + (-3) * x2 * sind(3 * phi) + 4 * nu * x2 * sind(3 * phi)) +
            (q3 + x3 + (-1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 +
            x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * (
            (-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * (5 + 4 * nu * ((
              - 3) + 2 * nu)) * y3p + (13 + 4 * nu * ((-7) + 4 * nu)) * x2 * cosd(phi) + ((-3) + 4 *
            nu) * x2 * cosd(3 * phi) + (-1) * (7 * q3 + 5 * x3 + (-20) * nu * (q3 + x3) + 16 *
            nu ^ 2 * (q3 + x3)) * sind(phi) + 2 * ((-3) + 4 * nu) * y2p * sind(2 * phi) + (3 *
            q3 + x3 + (-4) * nu * (q3 + x3)) * sind(3 * phi)) +(-4) * x3 * sind(phi) * (q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
            y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
            q3 + (-2) * y2p * cosd(phi) + 2 * (q3 + x3) * cosd(2 * phi) + (-1) * y2p * cosd(3 *
              phi) + y3p * sind(3 * phi)) +(-4) * x3 * (q3 + x3 + (-1) * y2p * cosd(phi) + y3p *
            sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 *
            y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) *
            sind(phi)) ^ (-2) * ((-2) * q3 * x2 * cosd(phi) + (3 * x2 * y2p + (q3 + (-1) * x3)
               * y3p) * cosd(2 * phi) + (-2) * y2p * y3p * cosd(3 * phi) + ((-1) * x2 * y2p + (
            q3 + x3) * y3p) * cosd(4 * phi) + ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (3 *
            q3 * y2p + x3 * y2p + (-1) * x2 * y3p) * sind(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3)
               ^ 2 + 2 * y2p ^ 2) * sind(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * sind(4 * phi)) +
            2 * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
            y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
               ^ (-1) * ((-2) * q3 * x2 * cosd(phi) + (3 * x2 * y2p + (q3 + (-1) * x3) * y3p) *
            cosd(2 * phi) + (-2) * y2p * y3p * cosd(3 * phi) + ((-1) * x2 * y2p + (q3 + x3) *
            y3p) * cosd(4 * phi) + ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * sind(phi) + (3 * q3 * y2p +
            x3 * y2p + (-1) * x2 * y3p) * sind(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 *
            y2p ^ 2) * sind(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * sind(4 * phi))) +xlogy((
              -1 / 2) * ((-2) + 4 * nu + cosd(2 * phi)) * sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) * q3 *
            x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p * cosd(
              phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 * y3p *
            sind(phi) + (-2) * x3 * y3p * sind(phi)) +xlogy((1 / 4) * ((-1) * (5 + (-20) *
              nu + 16 * nu ^ 2) * sind(phi) + (1 + (-4) * nu) * sind(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 *
            q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p *
            cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 *
            y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function I332d2(y2p::R, y3p::R) where R
            if phi ≈ 0 || phi ≈ 360
                ((1 / 2) * (8 * x3 * (x2 + (-1) * y3p) * (x2 ^ 2 + q3 * x3 + x3 ^ 2 + (-1) * x3 * y2p + (
                -2) * x2 * y3p + y3p ^ 2) * (q3 ^ 2 + x2 ^ 2 + x3 ^ 2 + 2 * q3 * (x3 + (-1) * y2p) + (
                -2) * x3 * y2p + y2p ^ 2 + (-2) * x2 * y3p + y3p ^ 2) ^ (-2) + (-16) * ((-1) + nu)
                   ^ 2 * (q3 + x3 + (-1) * y2p) * (x2 + (-1) * y3p) * (q3 ^ 2 + x2 ^ 2 + x3 ^ 2 + 2 * q3 *
                (x3 + (-1) * y2p) + (-2) * x3 * y2p + y2p ^ 2 + (-2) * x2 * y3p + y3p ^ 2) ^ (-1) +
                8 * ((-1) + nu) * (q3 + (-1) * x3 + (-1) * y2p) * (x2 + (-1) * y3p) * (q3 ^ 2 +
                x2 ^ 2 + x3 ^ 2 + 2 * x3 * y2p + y2p ^ 2 + (-2) * q3 * (x3 + y2p) + (-2) * x2 * y3p +
                y3p ^ 2) ^ (-1) + (-2) * ((-3) + 4 * nu) * (q3 + (-1) * x3 + (-1) * y2p) * (x2 + (
                -1) * y3p) * (q3 ^ 2 + x2 ^ 2 + x3 ^ 2 + 2 * x3 * y2p + y2p ^ 2 + (-2) * q3 * (x3 + y2p)
                +(-2) * x2 * y3p + y3p ^ 2) ^ (-1) + (10 * q3 + 6 * x3 + (-24) * nu * (q3 + x3) + 16 *
                nu ^ 2 * (q3 + x3) + (-2) * (5 + 4 * nu * ((-3) + 2 * nu)) * y2p) * (x2 + (-1) * y3p)
                   * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
                y2p + x2 * y3p)) ^ (-1) + 2 * x3 * ((-4) * x2 + 4 * y3p) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 *
                x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p)) ^ (-1)) +8 *
                ((-1) + nu) ^ 2 * atan((q3 + x3 + (-1) * y2p) * (x2 + (-1) * y3p) ^ (-1), 1) + (
                -4) * ((-1) + nu) * atan(((-1) * q3 + x3 + y2p) * ((-1) * x2 + y3p) ^ (-1), 1))
            elseif phi ≈ 180
                ((1 / 2) * (x2 + y3p) * ((-8) * ((-1) + nu) * (q3 + (-1) * x3 + y2p) * (q3 ^ 2 +
                x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + 2 * q3 * y2p + (-2) * x3 * y2p + y2p ^ 2 + 2 * x2 *
                y3p + y3p ^ 2) ^ (-1) + 2 * ((-3) + 4 * nu) * (q3 + (-1) * x3 + y2p) * (q3 ^ 2 +
                x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + 2 * q3 * y2p + (-2) * x3 * y2p + y2p ^ 2 + 2 * x2 *
                y3p + y3p ^ 2) ^ (-1) + (-8) * x3 * (x2 ^ 2 + q3 * x3 + x3 ^ 2 + x3 * y2p + 2 * x2 *
                y3p + y3p ^ 2) * (q3 ^ 2 + x2 ^ 2 + x3 ^ 2 + 2 * x3 * y2p + y2p ^ 2 + 2 * q3 * (x3 + y2p) +
                2 * x2 * y3p + y3p ^ 2) ^ (-2) + 8 * x3 * (q3 ^ 2 + x2 ^ 2 + x3 ^ 2 + 2 * x3 * y2p +
                y2p ^ 2 + 2 * q3 * (x3 + y2p) + 2 * x2 * y3p + y3p ^ 2) ^ (-1) + 16 * ((-1) + nu) ^ 2 *
                (q3 + x3 + y2p) * (q3 ^ 2 + x2 ^ 2 + x3 ^ 2 + 2 * x3 * y2p + y2p ^ 2 + 2 * q3 * (x3 + y2p) +
                2 * x2 * y3p + y3p ^ 2) ^ (-1) + (-2) * (5 * q3 + 3 * x3 + (-12) * nu * (q3 + x3) + 8 *
                nu ^ 2 * (q3 + x3) + (5 + 4 * nu * ((-3) + 2 * nu)) * y2p) * (q3 ^ 2 + x2 ^ 2 + x3 ^ 2 +
                2 * x3 * y2p + y2p ^ 2 + 2 * q3 * (x3 + y2p) + 2 * x2 * y3p + y3p ^ 2) ^ (-1)) +4 * ((
                  - 1) +nu) * atan((q3 + (-1) * x3 + y2p) * (x2 + y3p) ^ (-1), 1) + (-8) * ((-1) +
                nu) ^ 2 * atan((q3 + x3 + y2p) * (x2 + y3p) ^ (-1), 1))
            else
                (atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (y3p + (-1)
                   * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * cosd(phi) * (3 + (
                -4) * nu + cosd(2 * phi)) +(1 / 2) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (
                  -1) * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (
                -1), 1) * ((-1) * ((-13) + 4 * (7 + (-4) * nu) * nu) * cosd(phi) + (-1) * ((-3) +
                4 * nu) * cosd(3 * phi)) +(1 / 4) * ((-4) * (3 + (-4) * nu + cosd(2 * phi)) * ((-1)
                   * y3p + x2 * cosd(phi) + ((-1) * q3 + x3) * sind(phi)) * (q3 + (-1) * x3 + (-1) *
                y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 +
                y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2)
                   * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + 2 * (x2 + (-1) *
                y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 +
                x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(
                  phi) + (-2) * (x2 * y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6)
                   * y2p + 8 * nu * y2p + 2 * (q3 + (-1) * x3) * cosd(phi) * (2 + (-4) * nu + cosd(2 *
                    phi)) + 2 * x2 * (4 + (-4) * nu + cosd(2 * phi)) * sind(phi) + (-2) * y3p * sind(2 *
                      phi)) +2 * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)) * (q3 ^ 2 +
                x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 *
                y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * (
                (-2) * (5 + 4 * nu * ((-3) + 2 * nu)) * y2p + (7 * q3 + 5 * x3 + (-20) * nu * (q3 + x3)
                +16 * nu ^ 2 * (q3 + x3)) * cosd(phi) + (3 * q3 + x3 + (-4) * nu * (q3 + x3)) * cosd(
                  3 * phi) + (13 + 4 * nu * ((-7) + 4 * nu)) * x2 * sind(phi) + 2 * ((-3) + 4 * nu) *
                y3p * sind(2 * phi) + (3 + (-4) * nu) * x2 * sind(3 * phi)) +(-2) * (y3p + (-1) *
                x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (-2) * (q3 + x3 + (-1) * y2p * cosd(phi) +
                y3p * sind(phi)) * (1 + ((-1) * y2p + (q3 + x3) * cosd(phi) + x2 * sind(phi)) ^ 2 *
                (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (-2)) ^ (-1) * ((-2) * (
                5 + 4 * nu * ((-3) + 2 * nu)) * y3p + (13 + (-28) * nu + 16 * nu ^ 2) * x2 * cosd(phi)
                +2 * ((-3) + 4 * nu) * y3p * cosd(2 * phi) + (-1) * ((-3) + 4 * nu) * x2 * cosd(3 *
                  phi) + (-4) * nu * ((-5) + 4 * nu) * (q3 + x3) * sind(phi) + (-1) * (7 * q3 + 5 * x3)
                   * sind(phi) + (-1) * (3 * q3 + x3 + (-4) * nu * (q3 + x3)) * sind(3 * phi)) +4 *
                x3 * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
                y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
                   ^ (-1) * ((-2) * x2 * cosd(phi) + 3 * y3p * cosd(2 * phi) + (-2) * x2 * cosd(3 *
                  phi) + y3p * cosd(4 * phi) + (-2) * q3 * sind(phi) + y2p * sind(2 * phi) + y2p *
                sind(4 * phi)) +(-8) * x3 * (x2 + (-1) * y3p * cosd(phi) + (-1) * y2p * sind(phi)
                  ) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
                y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
                   ^ (-2) * ((-1) * ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * cosd(phi) + ((-1) * q3 * y2p +
                x3 * y2p + 3 * x2 * y3p) * cosd(2 * phi) + (-1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 *
                y3p ^ 2) * cosd(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * cosd(4 * phi) + (-2) *
                q3 * x2 * sind(phi) + (x2 * y2p + (3 * q3 + x3) * y3p) * sind(2 * phi) + (-2) *
                y2p * y3p * sind(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) * y3p) * sind(4 * phi))) +
                xlogy((1 / 2) * (4 + (-4) * nu + cosd(2 * phi)) * sind(phi), q3 ^ 2 + x2 ^ 2 + (-2) *
                  q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 * x3 * y2p *
                cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 *
                y3p * sind(phi) + (-2) * x3 * y3p * sind(phi)) + xlogy((1 / 4) * ((13 + 4 * nu * ((
                  - 7) + 4 * nu)) * sind(phi) + (3 + (-4) * nu) * sind(3 * phi)), q3 ^ 2 + x2 ^ 2 + 2 *
                q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) * x3 * y2p *
                cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) + 2 * q3 *
                y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
            end
        end

        function I332d3(y2p::R, y3p::R) where R
            (atan((y2p + ((-1) * q3 + x3) * cosd(phi) + (-1) * x2 * sind(phi)) * (y3p + (-1)
               * x2 * cosd(phi) + (q3 + (-1) * x3) * sind(phi)) ^ (-1), 1) * (3 + (-4) * nu + cosd(
              2 * phi)) * sind(phi) + (1 / 2) * atan((y2p + (-1) * (q3 + x3) * cosd(phi) + (-1)
               * x2 * sind(phi)) * (y3p + (-1) * x2 * cosd(phi) + (q3 + x3) * sind(phi)) ^ (-1),
            1) * ((-5) * sind(phi) + (-4) * nu * ((-5) + 4 * nu) * sind(phi) + (-1) * (1 + (-4)
               * nu) * sind(3 * phi)) +(1 / 2) * (2 * (3 + (-4) * nu + cosd(2 * phi)) * ((-1) *
            y3p + x2 * cosd(phi) + ((-1) * q3 + x3) * sind(phi)) * ((-1) * x2 + y3p * cosd(phi)
            +y2p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (
            -2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 * y2p + (-1)
               * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) + ((-1) * q3 + x3 + y2p * cosd(phi) + (
            -1) * y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + (-1) * x3 * y2p + x2 * y3p) * cosd(phi) + (-2) * (x2 *
            y2p + (-1) * q3 * y3p + x3 * y3p) * sind(phi)) ^ (-1) * ((-6) * y2p + 8 * nu *
            y2p + 2 * (q3 + (-1) * x3) * cosd(phi) * (2 + (-4) * nu + cosd(2 * phi)) + 2 * x2 * (
            4 + (-4) * nu + cosd(2 * phi)) * sind(phi) + (-2) * y3p * sind(2 * phi)) +(q3 + x3 + (
            -1) * y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 +
            y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) *
            x2 * y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-1) * ((-2) * (5 + 4 * nu * ((-3) + 2 *
            nu)) * y2p + (7 * q3 + 5 * x3 + (-20) * nu * (q3 + x3) + 16 * nu ^ 2 * (q3 + x3)) *
            cosd(phi) + (3 * q3 + x3 + (-4) * nu * (q3 + x3)) * cosd(3 * phi) + (13 + 4 * nu * ((
              - 7) + 4 * nu)) * x2 * sind(phi) + 2 * ((-3) + 4 * nu) * y3p * sind(2 * phi) + (3 + (
            -4) * nu) * x2 * sind(3 * phi)) +((-1) * x2 + y3p * cosd(phi) + y2p * sind(phi))
               * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (q3 * y2p + x3 *
            y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p) * sind(phi))
               ^ (-1) * (10 * y3p + (-24) * nu * y3p + 16 * nu ^ 2 * y3p + (-1) * (13 + (-28) *
            nu + 16 * nu ^ 2) * x2 * cosd(phi) + (-2) * ((-3) + 4 * nu) * y3p * cosd(2 * phi) + (
            -3) * x2 * cosd(3 * phi) + 4 * nu * x2 * cosd(3 * phi) + 7 * q3 * sind(phi) + (-20)
               * nu * q3 * sind(phi) + 16 * nu ^ 2 * q3 * sind(phi) + 5 * x3 * sind(phi) + (-20) *
            nu * x3 * sind(phi) + 16 * nu ^ 2 * x3 * sind(phi) + 3 * q3 * sind(3 * phi) + (-4) *
            nu * q3 * sind(3 * phi) + x3 * sind(3 * phi) + (-4) * nu * x3 * sind(3 * phi)) +4 *
            x3 * cosd(phi) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * (
            q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3) * y3p)
               * sind(phi)) ^ (-1) * (q3 + (-2) * (q3 + x3) * cosd(2 * phi) + y2p * cosd(3 * phi)
            +2 * y3p * sind(phi) + (-1) * y3p * sind(3 * phi)) +(-4) * x3 * (q3 + x3 + (-1) *
            y2p * cosd(phi) + y3p * sind(phi)) * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 +
            y3p ^ 2 + (-2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 *
            y2p + (q3 + x3) * y3p) * sind(phi)) ^ (-2) * ((-1) * ((-1) * q3 ^ 2 + x2 ^ 2 +
            x3 ^ 2) * cosd(phi) + ((-1) * q3 * y2p + x3 * y2p + 3 * x2 * y3p) * cosd(2 * phi) + (
            -1) * (x2 ^ 2 + (q3 + x3) ^ 2 + 2 * y3p ^ 2) * cosd(3 * phi) + ((q3 + x3) * y2p + x2 *
            y3p) * cosd(4 * phi) + (-2) * q3 * x2 * sind(phi) + (x2 * y2p + (3 * q3 + x3) * y3p)
               * sind(2 * phi) + (-2) * y2p * y3p * sind(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) *
            y3p) * sind(4 * phi)) +2 * (q3 ^ 2 + x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (
            -2) * (q3 * y2p + x3 * y2p + x2 * y3p) * cosd(phi) + 2 * ((-1) * x2 * y2p + (q3 + x3)
               * y3p) * sind(phi)) ^ (-1) * ((-1) * ((-1) * q3 ^ 2 + x2 ^ 2 + x3 ^ 2) * cosd(
              phi) + ((-1) * q3 * y2p + x3 * y2p + 3 * x2 * y3p) * cosd(2 * phi) + (-1) * (x2 ^ 2 +
            (q3 + x3) ^ 2 + 2 * y3p ^ 2) * cosd(3 * phi) + ((q3 + x3) * y2p + x2 * y3p) * cosd(4 *
              phi) + (-2) * q3 * x2 * sind(phi) + (x2 * y2p + (3 * q3 + x3) * y3p) * sind(2 * phi)
            +(-2) * y2p * y3p * sind(3 * phi) + (x2 * y2p + (-1) * (q3 + x3) * y3p) * sind(4 *
              phi))) +xlogy((-1 / 2) * cosd(phi) * (2 + (-4) * nu + cosd(2 * phi)), q3 ^ 2 +
                x2 ^ 2 + (-2) * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + 2 *
              x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) +
              2 * q3 * y3p * sind(phi) + (-2) * x3 * y3p * sind(phi)) + xlogy((1 / 4) * ((5 + (
                -20) * nu + 16 * nu ^ 2) * cosd(phi) + (1 + (-4) * nu) * cosd(3 * phi)), q3 ^ 2 +
            x2 ^ 2 + 2 * q3 * x3 + x3 ^ 2 + y2p ^ 2 + y3p ^ 2 + (-2) * q3 * y2p * cosd(phi) + (-2) *
            x3 * y2p * cosd(phi) + (-2) * x2 * y3p * cosd(phi) + (-2) * x2 * y2p * sind(phi) +
            2 * q3 * y3p * sind(phi) + 2 * x3 * y3p * sind(phi)))
        end

        function IU22(y2p::R, y3p::R) where R
            (1 / (8 * pi * G * (1 - nu)) * (
            sind(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I223d2(y2p, y3p)
            +2 * G * epsv23 * (I222d2(y2p, y3p) + I323d2(y2p, y3p))
            +(lambda * epsvkk + 2 * G * epsv33) * I322d2(y2p, y3p) )
            +cosd(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I222d2(y2p, y3p)
            +2 * G * epsv23 * (I322d2(y2p, y3p) - I223d2(y2p, y3p))
            -(lambda * epsvkk + 2 * G * epsv33) * I323d2(y2p, y3p))))
        end

        function IU23(y2p::R, y3p::R) where R
            (1 / (8 * pi * G * (1 - nu)) * (
            sind(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I223d3(y2p, y3p)
            +2 * G * epsv23 * (I222d3(y2p, y3p) + I323d3(y2p, y3p))
            +(lambda * epsvkk + 2 * G * epsv33) * I322d3(y2p, y3p) )
            +cosd(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I222d3(y2p, y3p)
            +2 * G * epsv23 * (I322d3(y2p, y3p) - I223d3(y2p, y3p))
            -(lambda * epsvkk + 2 * G * epsv33) * I323d3(y2p, y3p))))
        end

        function IU32(y2p::R, y3p::R) where R
            (1 / (8 * pi * G * (1 - nu)) * (
            sind(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I233d2(y2p, y3p)
            +2 * G * epsv23 * (I232d2(y2p, y3p) + I333d2(y2p, y3p))
            +(lambda * epsvkk + 2 * G * epsv33) * I332d2(y2p, y3p) )
            +cosd(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I232d2(y2p, y3p)
            +2 * G * epsv23 * (I332d2(y2p, y3p) - I233d2(y2p, y3p))
            -(lambda * epsvkk + 2 * G * epsv33) * I333d2(y2p, y3p) )))
        end

        function IU33(y2p::R, y3p::R) where R
            (1 / (8 * pi * G * (1 - nu)) * (
            sind(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I233d3(y2p, y3p)
            +2 * G * epsv23 * (I232d3(y2p, y3p) + I333d3(y2p, y3p))
            +(lambda * epsvkk + 2 * G * epsv33) * I332d3(y2p, y3p) )
            +cosd(phi) * ((lambda * epsvkk + 2 * G * epsv22) * I232d3(y2p, y3p)
            +2 * G * epsv23 * (I332d3(y2p, y3p) - I233d3(y2p, y3p))
            -(lambda * epsvkk + 2 * G * epsv33) * I333d3(y2p, y3p) ) ))
        end

        e22 = IU22(T / 2, W) - IU22(-T / 2, W) + IU22(-T / 2, zero(U)) - IU22(T / 2, zero(U))
        u23 = IU23(T / 2, W) - IU23(-T / 2, W) + IU23(-T / 2, zero(U)) - IU23(T / 2, zero(U))
        u32 = IU32(T / 2, W) - IU32(-T / 2, W) + IU32(-T / 2, zero(U)) - IU32(T / 2, zero(U))
        e33 = IU33(T / 2, W) - IU33(-T / 2, W) + IU33(-T / 2, zero(U)) - IU33(T / 2, zero(U))

        e22 = e22 - epsv22 * omega((x2 * sind(phi) - (x3 - q3) * cosd(phi)) / T) * S((x2 * cosd(phi) + (x3 - q3) * sind(phi)) / W)
        e23 = (u23 + u32) / 2 - epsv23 * omega((x2 * sind(phi) - (x3 - q3) * cosd(phi)) / T) * S((x2 * cosd(phi) + (x3 - q3) * sind(phi)) / W)
        e33 = e33 - epsv33 * omega((x2 * sind(phi) - (x3 - q3) * cosd(phi)) / T) * S((x2 * cosd(phi) + (x3 - q3) * sind(phi)) / W)

        ϵ[1] = e22
        ϵ[2] = e23
        ϵ[3] = e33
    end
end

function _stress_volinplane_rect!(σ::P, x2::U, x3::U, q2::U, q3::U, T::U, W::U, phi::U, epsv22p::U, epsv23p::U, epsv33p::U, G::U, nu::U) where {U, P<:AbstractVector{<:U}}
    _strain_volinplane_rect!(σ, x2, x3, q2, q3, T, W, phi, epsv22p, epsv23p, epsv33p, G, nu)
    lambda = G * 2 * nu / (1 - 2 * nu)
    e22, e23, e33 = σ
    s22 = lambda * (e22 + e33) + 2 * G * e22
    s23 = 2 * G * e23
    s33 = lambda * (e22 + e33) + 2 * G * e33
    σ[1], σ[2], σ[3] = s22, s23, s33
end

function _stress_volinplane_rect(x2::U, x3::U, q2::U, q3::U, T::U, W::U, phi::U, epsv22p::U, epsv23p::U, epsv33p::U, G::U, nu::U) where U
    σ = zeros(U, 3)
    _stress_volinplane_rect!(σ, x2, x3, q2, q3, T, W, phi, epsv22p, epsv23p, epsv33p, G, nu)
    return σ
end
