module RenewalDiD

using AutoHashEquals: @auto_hash_equals
using Compat: @compat
using DataFrames: DataFrame
using PrettyTables: pretty_table
using Random: AbstractRNG, default_rng
using StatsBase: Weights, mean, ordinalrank, quantile, sample
using Turing: @addlogprob!, @model
using Turing: Beta, Distribution, Exponential, LogNormal, Normal
using Turing: arraydist, cdf, filldist, product_distribution, truncated 

# re-export from `DataFrames`
export DataFrame
# interventionmatrix.jl
export InterventionMatrix
# generationinterval.jl
export g_covid, g_seir, generationtime, testgenerationtime, vectorg_seir
# simulations.jl
export packsimulations, packsimulationtuple, runsimulation, simulationcases, simulationu0
# fittingparameters.jl
export RenewalDiDData, RenewalDiDDataUnlimitedPopn, RenewalDiDPriors
export expectedseedcases, logsumexp, renewaldid
# processparameters.jl
export nunique, quantilerenewaldidinfections, rankvalues, samplerenewaldidinfections

# Plotting.jl 
@compat public plotmodel, plotmodel!
@compat public plotmodeldata, plotmodeldata!
@compat public plotmodelintervention, plotmodelintervention!
@compat public plotmodeloutput, plotmodeloutput!
@compat public traceplot, traceplot!, trplot, trplot!
@compat public tracerankplot, tracerankplot!

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            plottinghint = "\nHint: Plotting functions from `RenewalDiD` or \
                `RenewalDiD.Plotting` are only available after `using CairoMakie`"
            plottingfunctionslist = [
                plotmodel, plotmodel!,
                plotmodeldata, plotmodeldata!,
                plotmodelintervention, plotmodelintervention!,
                plotmodeloutput, plotmodeloutput!,
                traceplot, traceplot!, trplot, trplot!, 
                tracerankplot, tracerankplot!,
            ]

            if exc.f in plottingfunctionslist && 
                isnothing(Base.get_extension(RenewalDiD, :RenewalDiDCairoMakieExt))

                print(io, plottinghint)
            end
        end
    end
    return nothing
end

struct Automatic end  # not exported
const automatic = Automatic()  # not exported

# allow default values to be equal across multiple functions. In final version likely to 
# replace these constants with the values directly in the functions. 
const DEFAULT_SEEDMATRIX_HEIGHT = 7  # not exported 
const DEFAULT_SEEDMATRIX_MINVALUE = 0.5  # not exported

include("interventionmatrix.jl")
include("generationinterval.jl")
include("simulations.jl")
include("fittingparameters.jl")
include("processparameters.jl")

# submodules
include("Plotting.jl")
include("FittedParameterTestFunctions.jl")

end  # module RenewalDiD
