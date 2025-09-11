
using RenewalDiD
using Documenter
using Test

DocMeta.setdocmeta!(RenewalDiD, :DocTestSetup, :(using RenewalDiD, Turing); recursive=true)
doctest(RenewalDiD)
