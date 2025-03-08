using Test

#push!(LOAD_PATH, expanduser(".")) # Assumed to be ran from the repl folder.

#include("ArrayRadiation.jl")
using ArrayRadiation


include("doctest.jl")
include("test_gain_calculation.jl")