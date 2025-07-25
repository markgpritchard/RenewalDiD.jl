module RenewalDiD

using PrettyTables: pretty_table
using Random: AbstractRNG, default_rng
using Reexport: @reexport
using StatsBase: Weights, mean, ordinalrank, quantile, sample
using Turing: @model, Exponential, Normal, arraydist, filldist
using UnPack: @unpack

@reexport using DataFrames: DataFrame

struct Automatic end  # not exported
const automatic = Automatic()  # not exported

include("interventionmatrix.jl")
include("generationinterval.jl")
include("simulations.jl")
include("fittingparameters.jl")
include("processparameters.jl")

## interventionmatrix.jl
export InterventionMatrix
## generationinterval.jl
export g_covid, g_seir, generationtime, testgenerationtime, vectorg_seir
## simulations.jl
export packsimulations, packsimulationtuple, runsimulation, simulationcases, simulationu0
## fittingparameters.jl
export packdata, packpriors, renewaldid, renewaldid_tracksusceptibles
## processparameters.jl
export nunique, quantilerenewaldidinfections, rankvalues, samplerenewaldidinfections

## submodules
include("Plotting.jl")
include("FittedParameterTestFunctions.jl")

end  # module RenewalDiD
