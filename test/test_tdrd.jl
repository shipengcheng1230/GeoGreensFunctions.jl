using Test
using DelimitedFiles

@testset "Triangle Dislocation" begin
    P1 = [-1, -1, -5] * 1.0
    P2 = [1, -1, -5] * 1.0
    P3 = [-1, 1, -4] * 1.0
    ν = 0.25
    ss, ds, ts = 1.0, -1.0, 2.0

    xyz = [
        [-1/3, -1/3, -3.0],
        [-1/3, -1/3, -14/3],
        [-1/3, -1/3, -6.0],
        [7.0, -1.0, -5.0],
        [-7.0, -1.0, -5.],
        [-1.0, -3.0, -6.0],
        [-1.0, 3.0, -3.0],
        [3.0, -3.0, -6.0],
        [-3.0, 3.0, -3.0],
        [-1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
        [-1.0, -1.0, -8.0],
        [-1.0, 1.0, -8.0],
        [1.0, -1.0, -8.0],
        ]

    @testset "displacement" begin
        u = map(x -> disp_tri3_hs(x[1], x[2], x[3], P1, P2, P3, ss, ds, ts, ν)[1], xyz)
        # validation data from Table 3, Nikkhoo & Walter, 2015
        u_truth = [
            0.0352311877319734, -0.509465745232405, -0.0450664903903138,
            -0.00230579292792908, 0.00401472894583963, 0.00483219740842196,
            0.00261498816660580, -0.00498017723062124, 0.00469791958297159,
            0.00147786823004841, 0.0037484457568832, 0.0166558386729541,
            -0.00656347406353817, -0.0105680479573571, -0.00213929091658054,
            ]
        @test u ≈ u_truth
    end

    @testset "strain ϵₓₓ" begin
        ϵ = map(x -> strain_tri3_hs(x[1], x[2], x[3], P1, P2, P3, ss, ds, ts, 2.0, 2.0)[1], xyz)
        # validation data from Table 3, Nikkhoo & Walter, 2015
        ϵ_truth = [
            0.048104700525518, -0.244188978214975, 0.0546831404832553,
            0.000829157341339727, 0.00114439668841158, -0.00386292388925956,
            -0.00243788640223540, 0.000706397690338731, 0.000211254167350266,
            0.00650800501584133, 0.000922452413344460, 0.00441202690885827,
            0.00330232019558791, 0.00876398663844928, -0.000914111766849476,
            ]
        @test ϵ ≈ ϵ_truth
    end
end

@testset "consistency between rectangular and triangular" begin
    # settings are adopted from Figure 11, Nikkhoo & Walter, 2015
    l, w = 1.5, 0.75
    dip, strike = 30.0, 10.0
    depth = 2.0
    depth_obs = 5.0

    p1 = [l/2 * sind(strike), l/2 * cosd(strike), -depth]
    p2 = [-p1[1], -p1[2], -2.0]
    p4 = [p1[1] + w * cosd(dip) * cosd(strike), p1[2] - w * cosd(dip) * sind(strike), p1[3] - w * sind(dip)]
    p3 = [p4[1] - l * sind(strike), p4[2] - l * cosd(strike), p4[3]]
    p′ = map(deepcopy, [p1, p2, p3, p4])

    xs = range(-3.0, 3.0; step=1.0)
    ys = range(-3.0, 3.0; step=1.0)

    disl = [0.0, 1.0, 0.0]
    A = [cosd(90-strike) sind(90-strike) 0; -sind(90-strike) cosd(90-strike) 0; 0 0 1]
    A⁻¹ = inv(A)
    E = Matrix{Float64}(undef, 3, 3)

    for (x, y) in Iterators.product(xs, ys)
        obs = A * [x, y, -depth_obs]
        u_okada = dc3d(obs..., 2/3, depth, dip, [-l/2, l/2], [-w, 0.0], disl)
        ux, uy, uz = A⁻¹ * u_okada[1:3]
        u_td1 = disp_tri3_hs(x, y, -depth_obs, p1, p2, p3, disl..., 0.25)
        u_td2 = disp_tri3_hs(x, y, -depth_obs, p3, p4, p1, disl..., 0.25)
        ux_td, uy_td, uz_td = map(+, u_td1, u_td2)
        @test ux_td ≈ ux && uy_td ≈ uy && uz_td ≈ uz
        E[1,1] = u_okada[4]
        E[2,2] = u_okada[8]
        E[3,3] = u_okada[12]
        E[1,2] = (u_okada[5] + u_okada[7]) / 2
        E[1,3] = (u_okada[6] + u_okada[10]) / 2
        E[2,3] = (u_okada[9] + u_okada[11]) / 2
        E[2,1] = E[1,2]
        E[3,1] = E[1,3]
        E[3,2] = E[2,3]
        E′ = A' * E * A
        ϵ_td1 = strain_tri3_hs(x, y, -depth_obs, p1, p2, p3, disl..., 1.0, 1.0)
        ϵ_td2 = strain_tri3_hs(x, y, -depth_obs, p3, p4, p1, disl..., 1.0, 1.0)
        ϵ = collect(ϵ_td1) + collect(ϵ_td2)
        @test ϵ[1] ≈ E′[1,1] && ϵ[2] ≈ E′[2,2] && ϵ[3] ≈ E′[3,3] && ϵ[4] ≈ E′[1,2] && ϵ[5] ≈ E′[1,3] && ϵ[6] ≈ E′[2,3]
    end
    @test p′ == [p1, p2, p3, p4] # test if pts coordinates have been changed
end

@testset "Full space solution" begin
    xs = collect(-10.0: 0.5: 10.0)
    ys = collect(-8.0: 0.5: 8.0)
    zs = 4.0

    xyz = Iterators.product(zs, ys, xs) |> collect
    Y = [q[2] for q in xyz] |> vec
    X = [q[3] for q in xyz] |> vec
    Z = [q[1] for q in xyz] |> vec

    @testset "displacement" begin
        u, v, w = disp_rect_fs(X, Y, Z, 3.0, 1.0, 0.0, 2.0, 1.0, 10.0, 75.0, 30.0, 20.0, 1.0, 0.15, 0.25, Val(:pc))
        # validation data generated via orginal Matlab code
        truth = readdlm(joinpath(@__DIR__, "data/test_disp_rect_fs.dat"), ' ', Float64)
        @test truth[:,1] ≈ u && truth[:,2] ≈ v && truth[:,3] ≈ w
    end
    @testset "stress" begin
        stress = stress_rect_fs(X, Y, Z, 3.0, 1.0, 0.0, 2.0, 1.0, 10.0, 75.0, 30.0, 20.0, 1.0, 0.15, 0.66e11, 0.33e11, Val(:pc))
        # validation data generated via orginal Matlab code
        truth = readdlm(joinpath(@__DIR__, "data/test_stress_rect_fs.dat"), ' ', Float64)
        for i ∈ 1: 6
            @test truth[:,i] ≈ stress[i]
        end
    end
end
