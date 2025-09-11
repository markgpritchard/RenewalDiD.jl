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

"""
    RenewalDiD.plotmodelR0(A, [t; <keyword arguments>])

Plot changes in estimated R₀ over time.

Note that a `Makie.jl` backend must be loaded for this function to operate.

# Arguments
- `A`: `Array` of values. The third dimension is assume to represent different quantiles, 
    with the midpoint representing the median
- `t`: time (defaults to `axes(A, 1)`)

Keyword arguments are passed to the `Makie.jl` functions.

# Examples
```jldoctest
julia> import CairoMakie

julia> using StableRNGs

julia> rng = StableRNG(1000);

julia> data = RenewalDiD.testsimulation(rng);

julia> model = renewaldid(data, g_seir, RenewalDiDPriors(); mu=0.2, kappa=0.5);

julia> df = RenewalDiD.testdataframe(
       rng;
       nchains=2, niterations=5, ngroups=3, ntimes=10, nseeds=7,
       );
    
julia> A = samplerenewaldidinfections(rng, model, df);

julia> qs = quantilerenewaldidinfections(A.R0s, [0.05, 0.5, 0.95]);

julia> RenewalDiD.plotmodelR0(qs);
```
"""
function plotmodelR0 end 

"""
    RenewalDiD.plotmodelR0!(f, A, [t; <keyword arguments>])
    RenewalDiD.plotmodelR0!(axs, A, [t; <keyword arguments>])

Plot changes in estimated R₀ over time.

Mutating version of `RenewalDiD.plotmodelR0`. 

`f` can be a `Figure` or `GridLayout`, and `axs` can be an `Array` of `Axis`. See 
    `RenewalDiD.plotmodelR0` for details of other arguments.

# Examples
```jldoctest
julia> using CairoMakie, StableRNGs

julia> rng = StableRNG(1000);

julia> data = RenewalDiD.testsimulation(rng);

julia> model = renewaldid(data, g_seir, RenewalDiDPriors(); mu=0.2, kappa=0.5);

julia> df = RenewalDiD.testdataframe(rng; nchains=2, niterations=5, ngroups=3, ntimes=10, nseeds=7,);

julia> A = samplerenewaldidinfections(rng, model, df);

julia> qs = quantilerenewaldidinfections(A.R0s, [0.05, 0.5, 0.95]);

julia> fig1 = Figure();

julia> RenewalDiD.plotmodelR0!(fig1, qs);

julia> fig2 = Figure()

julia> axs = [Axis(fig2[1, i]) for i in 1:3];

julia> RenewalDiD.plotmodelR0!(axs, qs);
```
"""
function plotmodelR0! end

"""
    RenewalDiD.traceplot(df::DataFrame[, variable; <keyword arguments>])

Create a trace plot from Markov chain Monte Carlo samples.

Note that a `Makie.jl` backend must be loaded for this function to operate.

# Arguments
- `df`: DataFrame of samples 
- `variable`: String or Symbol or Array of these representing which columns to plot (default 
    is to plot all columns)

## Keyword arguments 
- `nplots`: number of columns to plot (not used with `variable` argument)
- `cols`: number of columns to arrange the plots into
Other keyword arguments are passed to the `Makie.jl` functions.

# Examples
```jldoctest
julia> import CairoMakie

julia> using StableRNGs

julia> rng = StableRNG(100);

julia> df = RenewalDiD.testdataframe(
       rng;
       nchains=4, niterations=2000, ngroups=2, ntimes=2, nseeds=5, 
       );

julia> RenewalDiD.traceplot(df, Symbol("tau[1]"));

julia> RenewalDiD.traceplot(df; ncols=5, nplots=25, size=(1000, 1000));
```
"""
function traceplot end

"""
    RenewalDiD.traceplot!(f, df::DataFrame[, variable; <keyword arguments>])
    RenewalDiD.traceplot!(axs, df::DataFrame[, variable; <keyword arguments>])

Create a trace plot from Markov chain Monte Carlo samples.

Mutating version of `RenewalDiD.traceplot`. 

`f` can be a `Figure` or `GridLayout`, and `axs` can be an `Array` of `Axis`. See 
    `RenewalDiD.traceplot` for details of other arguments.

# Examples
```jldoctest
julia> using CairoMakie, StableRNGs

julia> rng = StableRNG(100);

julia> df = RenewalDiD.testdataframe(
       rng;
       nchains=4, niterations=2000, ngroups=2, ntimes=2, nseeds=5, 
       );

julia> fig1 = Figure();

julia> RenewalDiD.traceplot!(fig1, df, Symbol("tau[1]"))

julia> fig2 = Figure();

julia> axs = [Axis(fig2[1, i]) for i in 1:3];

julia> RenewalDiD.traceplot!(axs, df, [:alpha, :sigma_gamma, Symbol("gammas_raw[1]")])
```
"""
function traceplot! end

"""
    RenewalDiD.tracerankplot(df::DataFrame[, variable; <keyword arguments>])

Create a ranked trace plot from Markov chain Monte Carlo samples.

Arguments are the same as for `RenewalDiD.traceplot`.

Note that a `Makie.jl` backend must be loaded for this function to operate.

# Examples
```jldoctest
julia> import CairoMakie

julia> using StableRNGs

julia> rng = StableRNG(100);

julia> df = RenewalDiD.testdataframe(
       rng;
       nchains=4, niterations=2000, ngroups=2, ntimes=2, nseeds=5, 
       );

julia> RenewalDiD.tracerankplot(df, Symbol("tau[1]"));

julia> RenewalDiD.tracerankplot(df; ncols=5, nplots=25, size=(1000, 1000));
```
"""
function tracerankplot end

"""
    RenewalDiD.tracerankplot!(f, df::DataFrame[, variable; <keyword arguments>])
    RenewalDiD.tracerankplot!(axs, df::DataFrame[, variable; <keyword arguments>])

Create a trace plot from Markov chain Monte Carlo samples.

Mutating version of `RenewalDiD.tracerankplot`. 

`f` can be a `Figure` or `GridLayout`, and `axs` can be an `Array` of `Axis`. See 
    `RenewalDiD.traceplot` for details of other arguments.

# Examples
```jldoctest
julia> using CairoMakie, StableRNGs

julia> rng = StableRNG(100);

julia> df = RenewalDiD.testdataframe(
       rng;
       nchains=4, niterations=2000, ngroups=2, ntimes=2, nseeds=5, 
       );

julia> fig1 = Figure();

julia> RenewalDiD.tracerankplot!(fig1, df, Symbol("tau[1]"))

julia> fig2 = Figure();

julia> axs = [Axis(fig2[1, i]) for i in 1:3];

julia> RenewalDiD.tracerankplot!(axs, df, [:alpha, :sigma_gamma, Symbol("gammas_raw[1]")])
```
"""
function tracerankplot! end

# `traceplot` is also exported by `Turing`
"""
    RenewalDiD.trplot(df::DataFrame[, variable; <keyword arguments>])

See `RenewalDiD.traceplot`. 
```
"""
trplot(args...; kwargs...) = traceplot(args...; kwargs...)

"""
    RenewalDiD.trplot!(f, df::DataFrame[, variable; <keyword arguments>])
    RenewalDiD.trplot!(axs, df::DataFrame[, variable; <keyword arguments>])

See `RenewalDiD.traceplot!`. 
```
"""
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
