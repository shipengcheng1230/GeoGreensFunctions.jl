module GeoGreensFunctions

const KERNELDIR = joinpath(@__DIR__, "funcs")
foreach(x -> include(joinpath(KERNELDIR, x)), filter!(x -> endswith(x, ".jl") && !startswith(x, "_"), readdir(KERNELDIR)))

end # module
