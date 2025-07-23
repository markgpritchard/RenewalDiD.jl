# plot outputs of RenewalDiD analysis

module RenewalDiDCairoMakie

using RenewalDiD
using CairoMakie: Axis, Figure, lines!

include("plottingfunctions/traceplot.jl")
include("plottingfunctions/tracerankplot.jl")


end  # module RenewalDiDCairoMakie
