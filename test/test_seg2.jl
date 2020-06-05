using Test

@testset "Antiplane line dislocation" begin
    u1 = dc3d(1.0, -1.0, -1.0, 2/3, 0.0, 90.0, [-1e6, 1e6], [-4.0, -3.0], [1.0, 0.0, 0.0])
    u2 = disp_antiplane(1.0, -1.0, 0.0, 4.0, 3.0, 1.0)
    @test u1[1] ≈ u2 # strike displacement
    G = 1.0
    σxy1 = G * (u1[5] + u1[7])
    σxz1 = G * (u1[6] + u1[10])
    σxy2, σxz2 = stress_antiplane(1.0, -1.0, 0.0, 4.0, 3.0, G, 1.0)
    @test σxy1 ≈ -σxy2 atol=1e-6
    @test σxz1 ≈ σxz2 atol=1e-6
end
