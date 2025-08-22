# functions for plotting 

# to avoid loading the `CairoMakie` package for users who are not planning to use it, 
# definitions of plotting functions are given in `RenewalDiDCairoMakieExt`. Plotting
# functions can be accessed by users either as `RenewalDiD.plotmodel` or with
# `using RenewalDiD.Plotting; plotmodel`

function plotmodel end
function plotmodel! end
function plotmodeldata end 
function plotmodeldata! end
function plotmodelintervention end
function plotmodelintervention! end
function plotmodeloutput end 
function plotmodeloutput! end
function plotmodelR0 end 
function plotmodelR0! end
function traceplot end
function traceplot! end
function tracerankplot end
function tracerankplot! end

# `traceplot` is also exported by `Turing`
trplot(args...; kwargs...) = traceplot(args...; kwargs...)
trplot!(args...; kwargs...) = traceplot!(args...; kwargs...)

module Plotting 

using RenewalDiD: plotmodel, plotmodel!
using RenewalDiD: plotmodeldata, plotmodeldata!
using RenewalDiD: plotmodelintervention, plotmodelintervention!
using RenewalDiD: plotmodeloutput, plotmodeloutput!
using RenewalDiD: plotmodelR0, plotmodelR0!
using RenewalDiD: traceplot, traceplot!, trplot, trplot!
using RenewalDiD: tracerankplot, tracerankplot!

export plotmodel, plotmodel!
export plotmodeldata, plotmodeldata!
export plotmodelintervention, plotmodelintervention!
export plotmodeloutput, plotmodeloutput!
export plotmodelR0, plotmodelR0!
export traceplot, traceplot!, trplot, trplot!
export tracerankplot, tracerankplot!

end  # module Plotting
