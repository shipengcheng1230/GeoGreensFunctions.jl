using Test
using DelimitedFiles
using Base.Iterators
using FastGaussQuadrature

# Those corresponding verifiable data are obtained from orginal matlab functions at
# https://bitbucket.org/sbarbot/bssa-2016237/src/master/
# https://bitbucket.org/sbarbot/bssa-2018058/src/default/
@testset "Volume Hex8" begin
    epsv11 = 11e-6
    epsv12 = 5e-6
    epsv13 = 6e-6
    epsv22 = 7e-6
    epsv23 = 8e-6
    epsv33 = 13e-6
    L = 40e3
    W = 20e3
    T = 10e3
    q1 = -10e3
    q2 = 10e3
    q3 = 1e3
    theta = 40.0
    G = 1.0
    nu = 0.25

    x1s = range(-10000.0, stop=10000.0, step=1000.0)
    x2s = range(-10000.0, stop=10000.0, step=1000.0)
    x3s = range(-30.0, stop=-20.0, step=5.0) # in real case should not be negative
    xxs = product(x1s, x2s, x3s)

    @testset "Displacement" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_disp_vol_hex8.dat"), ' ', Float64)
        for (i, x) in enumerate(xxs)
            u = _disp_vol_hex8(x..., q1, q2, q3, L, T, W, theta, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
            @test u ≈ u_truth[i,:]
        end
    end

    @testset "Stress" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_stress_vol_hex8.dat"), ' ', Float64)
        for (i, x) in enumerate(xxs)
            u = _stress_vol_hex8(x..., q1, q2, q3, L, T, W, theta, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
            @test ifelse(all(map(isnan, u)) && all(map(isnan, u_truth[i,:])), true, u ≈ u_truth[i,:])
        end
    end
end

@testset "Volume Tet4" begin
    epsv11 = 11e0
    epsv12 = 5e0
    epsv13 = 6e0
    epsv22 = 7e0
    epsv23 = 8e0
    epsv33 = 13e0
    G = 1.0
    nu = 0.25
    A = [1.0, 2.0, 3.0]
    B = [2.0, 5.0, 4.0]
    C = [11.0, 13.0, 12.0]
    D = [6.0, 5.0, 4.0]

    x1s = range(-100.0, stop=100.0, step=20.0)
    x2s = range(-100.0, stop=100.0, step=20.0)
    x3s = range(20.0, stop=30.0, step=5.0)
    xxs = product(x1s, x2s, x3s)
    qd = gausslegendre(15)

    @testset "Displacement" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_disp_vol_tet4.dat"), ' ', Float64)
        for (i, x) in enumerate(xxs)
            u = _disp_vol_tet4(qd, x..., A, B, C, D, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, nu)
            @test u ≈ u_truth[i, :]
        end
    end

    @testset "Stress" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_stress_vol_tet4.dat"), ' ', Float64)
        for (i, x) in enumerate(xxs)
            u = _stress_vol_tet4(qd, x..., A, B, C, D, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
            @test u ≈ u_truth[i, :]
        end
    end
end

@testset "Volume Inplane rect" begin
    epsv22 = 7e-6
    epsv23 = 8e-6
    epsv33 = 9e-6
    T = 10e3
    W = 20e3
    q2 = 15e3
    q3 = 1e-3
    phi = 0.0
    G = 1.0
    nu = 0.25

    x3 = range(-10000.0, stop=-1000.0, step=1000)
    x2 = range(-10000.0, stop=10000.0, step=2000)
    xxs = product(x2, x3)

    @testset "Stress" begin
        u_truth = readdlm(joinpath(@__DIR__, "data", "test_stress_volinplane_rect.dat"), ' ', Float64)
        for (i, x) in enumerate(xxs)
            u = _stress_volinplane_rect(x..., q2, q3, T, W, phi, epsv22, epsv23, epsv33, G, nu)
            @test u ≈ u_truth[i,:]
        end
    end
end

@testset "Coordinates: ENU <-> NED" begin
    x, y, z, qx, qy, qz, L, T, W, θ, ϵxx, ϵxy, ϵxz, ϵyy, ϵyz, ϵzz, G = rand(17)
    ν = 0.2
    u1 = disp_vol_hex8(x, y, -z, qx, qy, -qz, T, L, W, θ, ϵxx, ϵxy, ϵxz, ϵyy, ϵyz, ϵzz, G, ν)
    u2 = _disp_vol_hex8(y, x, z, qy, qx, qz, L, T, W, θ, ϵyy, ϵxy, -ϵyz, ϵxx, -ϵxz, ϵzz, G, ν)
    @test u1[1] ≈ u2[2] && u1[2] ≈ u2[1] && u1[3] ≈ -u2[3]

    u1 = stress_vol_hex8(x, y, -z, qx, qy, -qz, T, L, W, θ, ϵxx, ϵxy, ϵxz, ϵyy, ϵyz, ϵzz, G, ν)
    u2 = _stress_vol_hex8(y, x, z, qy, qx, qz, L, T, W, θ, ϵyy, ϵxy, -ϵyz, ϵxx, -ϵxz, ϵzz, G, ν)
    @test u1[1] ≈ u2[4] && u1[2] ≈ u2[2] && u1[3] ≈ -u2[5] && u1[4] ≈ u2[1] && u1[5] ≈ -u2[3] && u1[6] ≈ u2[6]

    u1 = stress_volinplane_rect(x, -z, qx, -qz, T, W, θ, ϵxx, ϵxz, ϵzz, G, ν)
    u2 = _stress_volinplane_rect(x, z, qx, qz, T, W, θ, ϵxx, -ϵxz, ϵzz, G, ν)
    @test u1[1] ≈ u2[1] && u1[2] ≈ -u2[2] && u1[3] ≈ u2[3]

    A, B, C, D = [rand(3) for _ in 1: 4]
    A2, B2, C2, D2 = map(x -> [x[2], x[1], -x[3]], [A, B, C, D])
    qd = gausslegendre(15)

    u1 = disp_vol_tet4(qd, x, y, -z, A2, B2, C2, D2, ϵxx, ϵxy, ϵxz, ϵyy, ϵyz, ϵzz, ν)
    u2 = _disp_vol_tet4(qd, y, x, z, A, B, C, D, ϵyy, ϵxy, -ϵyz, ϵxx, -ϵxz, ϵzz, ν)
    @test u1[1] ≈ u2[2] && u1[2] ≈ u2[1] && u1[3] ≈ -u2[3]

    u1 = stress_vol_tet4(qd, x, y, -z, A2, B2, C2, D2, ϵxx, ϵxy, ϵxz, ϵyy, ϵyz, ϵzz, G, ν)
    u2 = _stress_vol_tet4(qd, y, x, z, A, B, C, D, ϵyy, ϵxy, -ϵyz, ϵxx, -ϵxz, ϵzz, G, ν)
    @test u1[1] ≈ u2[4] && u1[2] ≈ u2[2] && u1[3] ≈ -u2[5] && u1[4] ≈ u2[1] && u1[5] ≈ -u2[3] && u1[6] ≈ u2[6]
end