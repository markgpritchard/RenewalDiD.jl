using Documenter
using RenewalDiD

DocMeta.setdocmeta!(RenewalDiD, :DocTestSetup, :(using RenewalDiD, Turing); recursive=true)
makedocs(sitename="RenewalDiD documentation")
