# plot outputs of RenewalDiD analysis

module RenewalDiDCairoMakieExt

using RenewalDiD
using RenewalDiD.Plotting
using CairoMakie: Axis, Figure, lines!

include("plottingfunctions/traceplot.jl")
include("plottingfunctions/tracerankplot.jl")

end  # module RenewalDiDCairoMakieExt
