module RenewalDiD

import NaNMath
import PrettyTables  # can remove this once support for version 2 no longer needed

using Compat: @compat
using DataFrames: DataFrame, insertcols!
using Distributions: Beta, Distribution, Exponential, LogNormal, Normal
using Distributions: cdf, product_distribution, truncated
using DistributionsAD: arraydist, filldist
using DynamicPPL: @addlogprob!, @model, Model
using MCMCChains: Chains
using Random: AbstractRNG, default_rng
using StatsBase: Weights, coef, coefnames, mean, ordinalrank, quantile, sample

@static if pkgversion(PrettyTables).major == 2
    using PrettyTables: pretty_table
else 
    using PrettyTables: @text__no_horizontal_lines, TextTableFormat, pretty_table
end

# re-export from `DataFrames`
export DataFrame
# interventionarrays.jl
export AbstractInterventionArray
export AbstractInterventionVector, InterventionVector
export AbstractInterventionMatrix, InterventionMatrix
export AbstractInterventionArray3, InterventionArray
export interventioncat, modifiedduration
# generationinterval.jl
export g_covid, g_seir, generationtime, testgenerationtime, vectorg_seir
# abstractrenewaldiddata.jl
export AbstractRenewalDiDData, RenewalDiDData, SimulationData
# simulations.jl
export packsimulations, packsimulationtuple, runsimulation, simulationcases, simulationu0
# fittingparameters.jl
export RenewalDiDPriors
export expectedseedcases, renewaldid, renewaldidpredmodel, renewaldidpredmodelnointervention
# processparameters.jl
export RenewalDiDModel, nunique, predictedR_0, quantilerenewaldidinfections, rankvalues
# functionsfortests.jl
@compat public testdataframe, testmodel, testsimulation
# Plotting.jl 
@compat public plotmodel, plotmodel!
@compat public plotmodeldata, plotmodeldata!
@compat public plotmodelintervention, plotmodelintervention!
@compat public plotmodeloutput, plotmodeloutput!
@compat public plotmodelR0, plotmodelR0!
@compat public traceplot, traceplot!, trplot, trplot!
@compat public tracerankplot, tracerankplot!

struct Automatic end  # not exported
const automatic = Automatic()  # not exported

include("interventionarrays.jl")
include("generationinterval.jl")
include("abstractrenewaldiddata.jl")
include("simulations.jl")
include("fittingparameters.jl")
include("processparameters.jl")
include("functionsfortests.jl")

# submodules
include("plotting/Plotting.jl")

const RENEWALDIDPLOTTINGFUNCTIONSLIST = [
    plotmodel, 
    plotmodel!,
    plotmodeldata, 
    plotmodeldata!,
    plotmodelintervention, 
    plotmodelintervention!,
    plotmodeloutput, 
    plotmodeloutput!,
    plotmodelR0,
    plotmodelR0!,
    traceplot, 
    traceplot!, 
    trplot, 
    trplot!, 
    tracerankplot, 
    tracerankplot!,
]

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            RDMakieExtisnothing = isnothing(Base.get_extension(RenewalDiD, :RenewalDiDMakieExt))
            if exc.f in RENEWALDIDPLOTTINGFUNCTIONSLIST && RDMakieExtisnothing
                print(
                    io, 
                    "\nHint: Plotting functions from `RenewalDiD` or `RenewalDiD.Plotting` \
                        are only available after `using CairoMakie`"
                )
            end
        end
    end
    return nothing
end

end  # module RenewalDiD
