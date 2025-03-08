using Documenter
using ArrayRadiation

# Run doctests for ArrayRadiation.jl

DocMeta.setdocmeta!(ArrayRadiation, :DocTestSetup, :(using ArrayRadiation); recursive=true)
Documenter.doctest(ArrayRadiation)