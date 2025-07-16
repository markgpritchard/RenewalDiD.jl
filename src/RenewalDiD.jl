module RenewalDiD

using PrettyTables: pretty_table
using Random: AbstractRNG, default_rng
using StatsBase: Weights, sample

include("automatic.jl")
include("interventionmatrix.jl")
include("generationinterval.jl")
include("simulations.jl")

## interventionmatrix.jl
export InterventionMatrix
## generationinterval.jl
export g_covid, g_seir, generationtime, vectorg_seir
## simulations.jl
export runsimulation, seirrates, simulateday!

end
