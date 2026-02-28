push!(LOAD_PATH, "../src/")

using Documenter

# Running `julia --project docs/make.jl` can be very slow locally.
# To speed it up during development, one can use make_local.jl instead.
# The code below checks wether its being called from make_local.jl or not.
const LOCAL = get(ENV, "LOCAL", "false") == "true"

if LOCAL
    include("../src/ArrayRadiation.jl")
    using .ArrayRadiation
else
    using ArrayRadiation
    ENV["GKSwstype"] = "100"
end

DocMeta.setdocmeta!(ArrayRadiation, :DocTestSetup, :(using ArrayRadiation); recursive=true)


makedocs(
    modules=[ArrayRadiation],
    format=Documenter.HTML(),
    sitename="ArrayRadiation.jl",
    pages=Any[
        "index.md",
        "api_reference.md",
        "Examples"=>Any[
            "Examples/element_radiation_pattern.md",
            "Examples/1D_array.md",
            "Examples/2D_array.md",
            "Examples/2D_array_3D_plot.md",
            "Examples/irregular_array.md",
            "Examples/linear_array_polar_gain_pattern.md",
            "Examples/window_function.md",
            "Examples/monopulse_pattern.md",
            "Examples/beam_steering.md",],
        "coordinate_system.md",
    ],
    doctest=true,
)

deploydocs(
    repo="github.com/ErikBuer/ArrayRadiation.jl.git",
    push_preview=true,
)