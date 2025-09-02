module RenewalDiD

using AutoHashEquals: @auto_hash_equals
using Compat: @compat
using DataFrames: DataFrame, insertcols!
using DynamicPPL: @addlogprob!, @model, Model
using PrettyTables: pretty_table
using Random: AbstractRNG, default_rng
using StatsBase: Weights, coef, coefnames, mean, ordinalrank, quantile, sample
using Turing: Beta, Chains, Distribution, Exponential, LogNormal, Normal
using Turing: arraydist, cdf, filldist, product_distribution, truncated
using Turing.Optimisation: ModeResult

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
export AbstractRenewalDiDData, RenewalDiDData, RenewalDiDDataUnlimitedPopn, SimulationData
# simulations.jl
export packsimulations, packsimulationtuple, runsimulation, simulationcases, simulationu0
# fittingparameters.jl
export RenewalDiDPriors, expectedseedcases, renewaldid
# processparameters.jl
export SampledOutput
export map_DataFrame
export nunique
export quantilerenewaldidinfections
export rankvalues
export samplerenewaldidinfections

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

# allow default values to be equal across multiple functions. In final version likely to 
# replace these constants with the values directly in the functions. 
const DEFAULT_SEEDMATRIX_HEIGHT = 7  # not exported 
const DEFAULT_SEEDMATRIX_MINVALUE = 0.5  # not exported

include("interventionarrays.jl")
include("generationinterval.jl")
include("abstractrenewaldiddata.jl")
include("simulations.jl")
include("fittingparameters.jl")
include("processparameters.jl")

# submodules
include("plotting/Plotting.jl")
include("fittedparametertestfunctions/FittedParameterTestFunctions.jl")

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
