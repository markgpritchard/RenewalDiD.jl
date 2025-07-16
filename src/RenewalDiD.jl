module RenewalDiD

using PrettyTables: pretty_table
using Random: AbstractRNG, default_rng
using Reexport: @reexport
using StatsBase: Weights, sample
using Turing: @model, arraydist, filldist
using UnPack: @unpack

@reexport using Distributions

include("automatic.jl")
include("interventionmatrix.jl")
include("generationinterval.jl")
include("simulations.jl")
include("fittingparameters.jl")

## interventionmatrix.jl
export InterventionMatrix
## generationinterval.jl
export g_covid, g_seir, generationtime, vectorg_seir
## simulations.jl
export packsimulations, runsimulation, simulationcases
## fittingparameters.jl
export packdata, packpriors, renewaldid

end
