using DelimitedFiles
using Test
using Base.Iterators

# Those corresponding verifiable data are obtained from orginal fortran subroutines at
# http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html
@testset "Test Okada's displacement" begin
    @testset "Normal computations" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_dc3d.dat"), ',', Float64)
        pos = product(-40.: -30., 20.: 30., 10.: 15.)
        fu = (x) -> dc3d(reverse(x)..., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        u_cal = map(fu, pos) |> vec
        cache = GeoGreensFunctions.dc3d_cache(Float64)
        for (i, x) in enumerate(pos)
            dc3d(reverse(x)..., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.], cache)
            @test u_truth[i,:] ≈ u_cal[i] ≈ cache[1]
        end
    end

    @testset "Negative depth" begin
        u = dc3d(10., 20., 30., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        @test isapprox(zeros(Float64, 12), u, rtol=1e-8)
    end

    @testset "Singularity" begin
        u_singular = zeros(Float64, 12)

        u = dc3d(-80., 0., -50., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d(-70., 0., -10., 2. /3, 10., 90., [-80., 120.], [-20., 20.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d(-80., 20., -50., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        u_truth = [
            -32.60576975, 59.39436703, 1.61816100, -0.43682428, 1.66503564, 0.27751175,
            0.64126234, -1.05836342, 0.03630938, -0.16974972, 0.61800758, 0.22842501]
        @test isapprox(u_truth, u, rtol=1e-8)


        u = dc3d(-10., 20., -50., 2. /3, 50., 90., [-80., 120.], [-30., 25.], [200., -150., 100.])
        u_truth = [
            -56.04511679, 52.38458483, 40.87786501, -0.10140811, -0.09099896, 0.02371607,
            1.60327214, -0.63369117, -1.13609685, 0.05783881, 0.88665943, 0.21946899,
        ]
        @test isapprox(u_truth, u, rtol=1e-8)

        u = dc3d(-10., 20., 0., 2. /3, 0., 0., [-80., 120.], [-30., 25.], [200., -150., 100.])
        # take a look: https://github.com/JuliaLang/julia/issues/23376
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d(-100., 20., 0., 2. /3, 0., 0., [-80., -90.], [20., 25.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d(-80., 20., 0., 2. /3, 0., 0., [-80., 100.], [25., 35.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)
    end
end
