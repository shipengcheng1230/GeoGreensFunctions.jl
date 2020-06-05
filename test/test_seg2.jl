using Test

@testset "Antiplane line dislocation" begin
    u1 = dc3d(0.0, -1.0, -1.0, 2/3, 0.0, 90.0, [-1e4, 1e4], [-4.0, -3.0], [1.0, 0.0, 0.0])
    u2 = disp_antiplane_seg2(1.0, -1.0, 0.0, 4.0, 3.0, 1.0)
    @test u1[1] ≈ u2 atol=1e-6 # strike displacement
    G = rand()*100
    σxy1 = G * (u1[5] + u1[7])
    σxz1 = G * (u1[6] + u1[10])
    σxy2, σxz2 = stress_antiplane_seg2(1.0, -1.0, 0.0, 4.0, 3.0, G, 1.0)
    @test σxy1 ≈ -σxy2 atol=1e-6
    @test σxz1 ≈ σxz2 atol=1e-6
end

@testset "Inplane line dislocation" begin
    θ = rand()*90
    G = rand()
    λ = rand()
    α = (λ + G) / (λ + 2G)
    ν = λ / 2 / (λ + G)
    uo = dc3d(0.0, -1.0, -1.0, α, 0.0, θ, [-1e5, 1e5], [-4.0, -3.0], [0.0, 1.0, 0.0])
    u1, u2 = disp_inplane_seg2(1.0, -1.0, 4.0*cosd(θ), -4.0*sind(θ), 3.0*cosd(θ), -3.0*sind(θ), -cosd(θ), sind(θ), ν)
    @test uo[2] ≈ -u1 atol=1e-6
    @test uo[3] ≈ u2 atol=1e-6
    σ11, σ22, σ12 = stress_inplane_seg2(1.0, -1.0, 4.0*cosd(θ), -4.0*sind(θ), 3.0*cosd(θ), -3.0*sind(θ), -cosd(θ), sind(θ), G, ν)
    ukk = uo[4] + uo[8] + uo[12]
    σyy = λ * ukk + 2G * uo[8]
    σzz = λ * ukk + 2G * uo[12]
    σyz = G * (uo[9] + uo[11])
    @test σ11 ≈ σyy atol=1e-6
    @test σ22 ≈ σzz atol=1e-6
    @test σ12 ≈ -σyz atol=1e-6
end
