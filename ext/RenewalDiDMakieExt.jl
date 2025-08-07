# plot outputs of RenewalDiD analysis

module RenewalDiDMakieExt

using RenewalDiD
using RenewalDiD: Automatic, automatic
using Makie: Axis, Cycled, Figure, GridLayout, Label
using Makie: band!, lines!, linkaxes!, scatter!, vlines!

const FigOrGridLayout = Union{Figure, GridLayout} 
const StringOrSymbol = Union{<:AbstractString, Symbol} 
const StringOrSymbolOrMissing = Union{<:AbstractString, Symbol, Missing} 

include("makieext/keywords.jl")
include("makieext/plotnames.jl")
include("makieext/traceplot.jl")
include("makieext/tracerankplot.jl")
include("makieext/plotmodel.jl")

end  # module RenewalDiDMakieExt
