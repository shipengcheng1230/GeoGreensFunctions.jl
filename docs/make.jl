# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter
using GeoGreensFunctions

makedocs(
    doctest=false,
    modules = [GeoGreensFunctions],
    sitename = "GeoGreensFunctions",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "APIs" => "APIs.md",
    ],
)

deploydocs(
  repo = "github.com/shipengcheng1230/GeoGreensFunctions.git",
  target = "build",
)
