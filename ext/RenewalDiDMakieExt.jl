# plot outputs of RenewalDiD analysis

module RenewalDiDMakieExt

using RenewalDiD
using RenewalDiD.Plotting
using RenewalDiD: Automatic, automatic
using Makie: Axis, Cycled, Figure, GridLayout, band!, lines!, linkaxes!, scatter!, vlines! 

const FigOrGridLayout = Union{Figure, GridLayout}  # not exported

include("makieext/keywords.jl")
include("makieext/traceplot.jl")
include("makieext/tracerankplot.jl")
include("makieext/plotmodel.jl")

end  # module RenewalDiDMakieExt
