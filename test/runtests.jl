using GeoGreensFunctions

const TESTDIR = @__DIR__
const TESTFILES = filter(x -> startswith(x, "test_") && endswith(x, ".jl"), readdir(TESTDIR))
foreach(x -> include(joinpath(TESTDIR, x)), TESTFILES)
