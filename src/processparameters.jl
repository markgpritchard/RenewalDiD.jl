# run simulation with fitted parameters 

@auto_hash_equals struct SampledOutput{T, N} 
    output::Array{T, N}
    R0s::Array{T, N}

    function SampledOutput(output::Array{T, N}, R0s::Array{T, N}) where {T, N}
        return new{T, N}(output, R0s)
    end
end

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## number unique

"""
    nunique(vec)

Number of unique elements in a vector.

# Examples
```jldoctest
julia> nunique([1, 3, 2, 3])
3

julia> nunique(["a", "d", "d", "d"])
2
```
"""
nunique(vec) = length(unique(vec))

## data for tracerankplot 

"""
    rankvalues(df::DataFrame, variable; binsize=10)

Return rank of values from different chains for use in a trace rank plot.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(100);

julia> df = DataFrame(:chain => repeat([1, 2, 3]; inner=5), :iteration => repeat(1:5; \
    outer=3), :a => rand(rng, 15));  

julia> rankvalues(df, :a; binsize=3)
15-element Vector{Float64}:
 1.0
 1.0
 1.0
 2.0
 2.0
 2.6666666666666665
 2.6666666666666665
 2.6666666666666665
 2.0
 2.0
 2.3333333333333335
 2.3333333333333335
 2.3333333333333335
 2.0
 2.0
```
"""
function rankvalues(df::DataFrame, variable; binsize=10)
    vals = getproperty(df, variable)
    nchains = nunique(getproperty(df, :chain))
    niterations = round(Int, size(df, 1) / nchains)
    _rankvalueassertions(df, niterations)
    ranks = _unbinnedrankvalues(vals, nchains, niterations)
    binnedranks = _rankvalues(ranks, nchains, niterations, binsize)
    return binnedranks
end

function _unbinnedrankvalues(vals, nchains, niterations)
    ranks = zeros(Int, nchains * niterations)

    for i in 1:niterations
        inds = _rankvalueinds(nchains, niterations, i)
        ranks[inds] .= ordinalrank(vals[inds])
    end

    return ranks
end

function _rankvalues(ranks, nchains, niterations, binsize)
    binnedranks = zeros(nchains * niterations)
    
    bininds = _rankvaluebininds(nchains, niterations, binsize)
    for inds in bininds 
        binnedranks[inds] .= mean(ranks[inds])
    end

    return binnedranks
end

_rankvalueinds(nchains, niterations, i) = [(ch - 1) * niterations + i for ch in 1:nchains]

function _rankvaluebininds(nchains, niterations, binsize)
    inds = Vector{Int64}[]
    start = 1
    iterationend = niterations
    while start <= nchains * niterations
        ind = __rankvaluebininds(start, binsize, iterationend)
        push!(inds, ind)
        start = last(ind) + 1
        if last(ind) == iterationend 
            iterationend += niterations 
        end
    end
    return inds
end

function __rankvaluebininds(start, binsize, iterationend)
    stop = min(start + binsize - 1, iterationend)
    return start:stop
end

## versions of gamma and theta vector functions taking a `DataFrame` of fitted parameters

function _gammavec(df::DataFrame, i, ngroups)
    return _gammavec(  # call version in `fittingparameters.jl`
        [getproperty(df, Symbol("gammas_raw[$j]"))[i] for j in 1:(ngroups - 1)], 
        df.sigma_gamma[i]
    )
end 

function _thetavec(df::DataFrame, i, ntimes)
    return _thetavec(  # call version in `fittingparameters.jl`
        [getproperty(df, Symbol("thetas_raw[$j]"))[i] for j in 1:(ntimes - 1)], 
        df.sigma_theta[i]
    )
end  

function _tauvec(df, i, ninterventions)
    return [getproperty(df, Symbol("tau[$k]"))[i] for k in 1:ninterventions]
end  

function _mxmatrix(df, i, ngroups, ntimes, n_seeds)
    totaltimes = ntimes + n_seeds
    return [getproperty(df, Symbol("M_x[$t, $j]"))[i] for t in 1:totaltimes, j in 1:ngroups]
end

## take samples from DataFrame of fitted parameters and generate expected outcomes

"""
    samplerenewaldidinfections(g, df, data, indexes; <keyword arguments>)

Generate expected outcomes from samples in a `DataFrame`.

# Arguments 
- `g`: generation interval.
- `df::DataFrame`: outputs from parameter fitting.
- `data::AbstractRenewalDiDData`: data used for parameter fitting.
- `indexes::Union{<:AbstractVector{<:Integer}, <:Integer}=axes(df, 1)`: rows of DataFrame to 
    be used.
- Keyword arguments get passed to the generation interval function 


# Examples
```jldoctest
julia> using RenewalDiD.FittedParameterTestFunctions, StableRNGs

julia> rng = StableRNG(1000);

julia> data = testsimulation(rng);

julia> df = testdataframe(rng; nchains=2, niterations=5, ngroups=3, ntimes=10, nseeds=7);

julia> samplerenewaldidinfections(g_seir, df, data, 1; mu=0.2, kappa=0.5)
11×3 Matrix{Float64}:
 0.0        0.0       0.0
 0.114199   0.114506  0.029641
 0.136986   0.137571  0.0285553
 0.186847   0.188042  0.0288459
 0.266959   0.477386  0.0299742
 0.389018   0.828994  0.0316884
 0.573569   1.50836   0.0600483
 0.852888   2.78382   0.069359
 1.27559    5.15169   0.0831346
 1.91162    9.42535   0.101564
 2.85361   16.5956    0.125572
```
"""
function samplerenewaldidinfections(
    model, fittedparameters, indexes=automatic; 
    kwargs...
)
    return samplerenewaldidinfections(
        default_rng(), model, fittedparameters, indexes; 
        kwargs...
    )
end

function samplerenewaldidinfections(
    rng::AbstractRNG, 
    model, 
    fittedparameters::Tuple{DataFrame, <:Chains}, 
    indexes=automatic; 
    kwargs...
)
    return samplerenewaldidinfections(rng, model, fittedparameters[1], indexes; kwargs...)
end

function samplerenewaldidinfections(
    rng::AbstractRNG, model, fittedparameters::Chains, indexes=automatic; 
    kwargs...
)
    return samplerenewaldidinfections(
        rng, model, DataFrame(fittedparameters), indexes; 
        kwargs...
    )
end

function samplerenewaldidinfections(
    rng::AbstractRNG, model, fittedparameters::DataFrame, indexes=automatic; 
    repeatsamples=nothing, kwargs...
)
    return _samplerenewaldidinfections(
        rng, model, fittedparameters, indexes, repeatsamples; 
        kwargs...
    )
end

function samplerenewaldidinfections(
    ::AbstractRNG, ::Any, fittedparameters, x=nothing; 
    kwargs...
) 
    throw(_dftypeerror(fittedparameters))
    return nothing
end

function _samplerenewaldidinfections(
    rng, model, fittedparameters, ::Automatic, repeatsamples; 
    kwargs...
) 
    return _samplerenewaldidinfections(
        rng, model, fittedparameters, axes(fittedparameters, 1), repeatsamples; 
        kwargs...
    )
end

function _samplerenewaldidinfections(
    rng, model, fittedparameters, indexes, repeatsamples; 
    kwargs...
)
    model.f == _renewaldid || throw(ArgumentError("to do"))
    ngroups = _ngroups(model.args.interventions)
    ntimes = _ntimes(model.args.interventions)
    R0s = _samplerenewaldidinitialoutput(ntimes - 1, ngroups, indexes, repeatsamples)
    output = _samplerenewaldidinitialoutput(ntimes, ngroups, indexes, repeatsamples)
    _samplerenewaldidinfections!(
        rng, output, R0s, model, fittedparameters, indexes, repeatsamples; 
        kwargs...
    )
    return SampledOutput(output, R0s)
end

function _samplerenewaldidinitialoutput(ntimes, ngroups, ::Integer, ::Nothing)
    return zeros(ntimes + 1, ngroups)
end

function _samplerenewaldidinitialoutput(
    ntimes, ngroups, indexes::AbstractVector{<:Integer}, ::Nothing
)
    return zeros(ntimes + 1, ngroups, length(indexes))
end

function _samplerenewaldidinitialoutput(ntimes, ngroups, ::Integer, repeatsamples::Integer)
    return zeros(ntimes + 1, ngroups, repeatsamples)
end

function _samplerenewaldidinitialoutput(
    ntimes, ngroups, indexes::AbstractVector{<:Integer}, repeatsamples::Integer
)
    return zeros(ntimes + 1, ngroups, length(indexes) * repeatsamples)
end

function _samplerenewaldidinfections!(
    rng, output, R0s, model, fittedparameters, indexes::AbstractVector{<:Integer}, ::Nothing; 
    kwargs...
)
    for (r, i) in enumerate(indexes)
        _samplerenewaldidinfections!(
            rng, 
            (@view output[:, :, r]), 
            (@view R0s[:, :, r]), 
            model, 
            fittedparameters, 
            i, 
            nothing; 
            kwargs...
        )
    end
    return nothing
end

function _samplerenewaldidinfections!(
    rng,
    output, 
    R0s,
    model, 
    fittedparameters, 
    indexes::AbstractVector{<:Integer}, 
    repeatsamples::Number; 
    kwargs...
)
    r = 1
    for i in indexes
        for _ in 1:repeatsamples 
            _samplerenewaldidinfections!(
                rng, 
                (@view output[:, :, r]), 
                (@view R0s[:, :, r]), 
                model, 
                fittedparameters, 
                i, 
                nothing; 
                kwargs...
            )
            r += 1
        end
    end
    return nothing
end

function _samplerenewaldidinfections!(
    rng, output, R0s, model, fittedparameters, i::Integer, repeatsamples::Number; 
    kwargs...
)
    for r in 1:repeatsamples 
        _samplerenewaldidinfections!(
            rng, 
            (@view output[:, :, r]), 
            (@view R0s[:, :, r]), 
            model, 
            fittedparameters, 
            i, 
            nothing; 
            kwargs...
        )
    end
    return nothing
end

function _samplerenewaldidinfections!(
    rng, output, R0s, model, fittedparameters, i::Integer, ::Nothing; 
    kwargs...
)
    ngroups = _ngroups(model.args.interventions)
    ntimes = _ntimes(model.args.interventions)
    _samplerenewaldidinfectionsassertions(
        fittedparameters, model.args.expectedseedcases, i, ngroups, ntimes
    )
    alpha = fittedparameters.alpha[i]
    gammavec = _gammavec(fittedparameters, i, ngroups)
    thetavec = _thetavec(fittedparameters, i, ntimes)
    tau = _tauvec(fittedparameters, i, _ninterventions(model.args.interventions))
    psi = fittedparameters.psi[i]
    M_x = _mxmatrix(fittedparameters, i, ngroups, ntimes, model.args.n_seeds)
    R0s .= _predictedlogR_0(alpha, gammavec, thetavec, tau, model.args.interventions)
    T = Complex{typeof(R0s[1, 1])}
    predictedinfections = _infectionsmatrix(T, R0s, model.args.n_seeds)
    _infections!(
        model.args.g, 
        predictedinfections, 
        M_x, 
        R0s, 
        model.args.expectedseedcases, 
        model.args.Ns, 
        model.args.n_seeds; 
        model.defaults...
    )
    delayedinfections = _delayedinfections(
        T, predictedinfections, model.args.delaydistn, ngroups, ntimes, model.args.n_seeds
    )
    np = real.(delayedinfections[model.args.n_seeds:(model.args.n_seeds + ntimes), :]) .* psi
    isnan(maximum(np)) && return _halt_samplerenewaldidinfections(i)  # exit early 
    output .= _predictobservedinfections(rng, np, psi)
    return nothing
end

function _predictobservedinfections(rng, np::AbstractMatrix, psi)
    return [
        __predictobservedinfections(rng, np[t, g], psi) 
        for t in axes(np, 1), g in axes(np, 2)
    ]
end

function __predictobservedinfections(rng, np::Number, psi)
    v = rand(rng, Normal(np, sqrt(np * (1 - psi))))
    v > 0 || return zero(v)
    return v
end

function _halt_samplerenewaldidinfections(i) 
    @warn "Parameters in row $i give `NaN` predicted infections"
    return nothing 
end

## quantiles of sampled outputs

"""
    quantilerenewaldidinfections(A, q)

Calculate quantiles from an array of simulated outputs

# Arguments 
- `A`: array of outputs.
- `q`: vector of quantiles.

It is expected that an odd number of quantiles will be supplied, symmetrical around the 
    median (`0.5`)

# Examples
```jldoctest
julia> using RenewalDiD.FittedParameterTestFunctions, StableRNGs

julia> rng = StableRNG(1000);

julia> data = testsimulation(rng);

julia> df = testdataframe(rng; nchains=2, niterations=5, ngroups=3, ntimes=10, nseeds=7);

julia> A = samplerenewaldidinfections(g_seir, df, data; gamma=0.2, sigma=0.5);

julia> quantilerenewaldidinfections(A, [0.05, 0.5, 0.95])
11×3×3 Array{Float64, 3}:
[:, :, 1] =
 0.0        0.0        0.0
 0.0330109  0.0331442  0.0153384
 0.0424594  0.0426307  0.0159759
 0.0506862  0.0509582  0.0164732
 0.0640469  0.114019   0.0174583
 0.0843459  0.159901   0.018906
 0.113446   0.23254    0.0412647
 0.151611   0.334978   0.0491071
 0.0768479  0.108033   0.0195905
 0.100981   0.1344     0.0274011
 0.0        0.0        0.0397927

[:, :, 2] =
 0.0        0.0        0.0
 0.0979914  0.0982897  0.0616108
 0.119181   0.119662   0.0697299
 0.165515   0.166476   0.0985363
 0.328753   0.576721   0.177187
 0.535071   0.949926   0.38878
 0.839773   1.63126    0.824038
 1.4096     2.84745    1.44794
 1.70084    5.003      1.09146
 2.51385    8.74408    1.62129
 1.8444     9.59518    1.47483

[:, :, 3] =
  0.0        0.0        0.0
  0.289015   0.290399   0.202453
  0.401929   0.405079   0.25646
  0.660351   0.668552   0.371832
  1.18571    1.61497    0.579842
  2.27498    3.94403    0.950158
  6.11523   19.8733     2.91164
 13.4612    19.676      5.75056
 20.3066    47.5492     8.25488
 32.011     75.8905    14.145
 20.3503    47.2522    17.3482

julia> quantilerenewaldidinfections(A, [0.05, 0.5, 0.90]);
┌ Warning: (0.05, 0.9): other functions expect that credible intervals are symmetrical
└ @ RenewalDiD 

julia> quantilerenewaldidinfections(A, [0.05, 0.90]);
┌ Warning: [0.05, 0.9]: other functions expect an odd number of quantiles
└ @ RenewalDiD 
```
"""
function quantilerenewaldidinfections(SO::SampledOutput, q; mutewarnings=nothing)
    return SampledOutput(
        quantilerenewaldidinfections(SO.output, q; mutewarnings),
        quantilerenewaldidinfections(SO.R0s, q; mutewarnings)
    )
end

function quantilerenewaldidinfections(A::AbstractArray, q; mutewarnings=nothing)
    _quantilerenewaldidinfectionswarningset(A, q, mutewarnings)
    return _quantilerenewaldidinfections(A, q)
end

_quantilerenewaldidinfections(M::AbstractMatrix, ::Any) = M

function _quantilerenewaldidinfections(A::Array{T, 3}, q) where T
    outputnumbertype = typeof(zero(T) / 1)
    return __quantilerenewaldidinfections(outputnumbertype, A, q)
end

function __quantilerenewaldidinfections(outputnumbertype, A, q::Real)
    output = zeros(outputnumbertype, size(A, 1), size(A, 2))
    for t in axes(A, 1), j in axes(A, 2)
        output[t, j] = quantile((@view A[t, j, :]), q)
    end
    return output
end

function __quantilerenewaldidinfections(outputnumbertype, A, q::AbstractVector{<:Real}) 
    output = zeros(outputnumbertype, size(A, 1), size(A, 2), length(q))
    for (i, qi) in enumerate(q), t in axes(A, 1), j in axes(A, 2)
        output[t, j, i] = quantile((@view A[t, j, :]), qi)
    end
    return output
end

## list the intervention times

# other versions of this function are in `interventionmatrix.jl`
_interventionstarttimes(M::AbstractMatrix, i) = findfirst(x -> x == 1, M[:, i])

## DataFrame of mode estimates

function map_DataFrame(result::ModeResult)
    df = DataFrame()
    for (n, v) in zip(coefnames(result), coef(result))
        insertcols!(df, n => v)
    end
    return df
end


# Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Warning sets 

function _quantilerenewaldidinfectionswarningset(A, q, mutewarnings::Bool)
    if mutewarnings 
        return nothing 
    else 
        return _quantilerenewaldidinfectionswarningset(A, q, nothing)
    end
end

function _quantilerenewaldidinfectionswarningset(::AbstractMatrix, ::Any, ::Nothing)
    return _singlesamplequantilewarning()
end

function _quantilerenewaldidinfectionswarningset(A::AbstractArray{<:Any, 3}, q, ::Nothing) 
    size(A, 3) > 1 || _singlesamplequantilewarning()
    return __quantilerenewaldidinfectionswarningset(q)
end

__quantilerenewaldidinfectionswarningset(::Real) = nothing

function __quantilerenewaldidinfectionswarningset(q::AbstractVector)
    isodd(length(q)) || return _evenquantilevectorlengthwarning(q)
    _midpoint = convert(Int, length(q) / 2 + 0.5)
    q[_midpoint] == 0.5 || _medianmidpointwarning(q)
    qp = 0
    for (i, qi) in enumerate(q) 
        qi > qp || _descendingquantilewarning(qi, qp)
        qp = qi 
        i >= _midpoint && continue 
        qo = q[(length(q) + 1 - i)]
        qi ≈ 1 - qo || _assymetricquantileswarning(qi, qo)
    end
    return nothing
end

## Warning messages 

function _assymetricquantileswarning(qi, qo)
    @warn "($qi, $qo): other functions expect that credible intervals are symmetrical"
    return nothing
end

function _descendingquantilewarning(qi, qp)
    @warn "$qi < $qp: other functions expect quantiles are calculated in ascending order"
    return nothing
end

function _evenquantilevectorlengthwarning(q)
    @warn "$q: other functions expect an odd number of quantiles"
    return nothing 
end

function _medianmidpointwarning(q)
    @warn "$q: other functions expect the middle quantile is the median"
    return nothing
end

function _singlesamplequantilewarning()
    wm = "only a single sample provided to `quantilerenewaldidinfections`; input returned \
        unchanged for any value of `q`"
    @warn wm
    return nothing
end


# Assertions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _rankvalueassertions(df, niterations)
    chainnumbers = unique(df.chain)
    df.chain == repeat(chainnumbers; inner = niterations) || throw(_rankvaluechainerror())
    return nothing
end

function _samplerenewaldidinfectionsassertions(df, seedmatrix, i, ngroups, ntimes)
    _ngroups(seedmatrix) == ngroups || throw(_seedwidthdimensionerror())
    M_xmaxt = _ntimes(seedmatrix) + ntimes
    "M_x[$M_xmaxt, 1]" in names(df) || throw(_dfMxdimensiontoosmallerror())
    "M_x[$(M_xmaxt + 1), 1]" ∉ names(df) || throw(_dfMxdimensiontoolargererror())
    return _samplerenewaldidinfectionsindexassertions(df, i)
end

function _samplerenewaldidinfectionsindexassertions(df, i::Integer)
    i >= 1 || throw(BoundsError(df, [i, 1]))
    i <= size(df, 1) || throw(BoundsError(df, [i, 1]))
    return nothing
end

function _samplerenewaldidinfectionsindexassertions(df, indexes::AbstractVector{<:Integer})
    minimum(indexes) >= 1 || throw(BoundsError(df, [indexes, 1]))
    maximum(indexes) <= size(df, 1) || throw(BoundsError(df, [indexes, 1]))
    return nothing
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _dfMxdimensiontoosmallerror()
    m = "`df` must contain values of M_x for all times up to `_ntimes(seedmatrix) + ntimes`"
    return DimensionMismatch(m)
end

function _dfMxdimensiontoolargererror()
    m = "`df` must not contain values M_x for times beyond `_ntimes(seedmatrix) + ntimes`"
    return DimensionMismatch(m)
end

function _dftypeerror(::T) where T 
    m = "`fittedparameters` must be a `DataFrame`, `Chains` or `Tuple{DataFrame, <:Chains}`"
    return ArgumentError("fittedparameters::$T: " * m)
end

_rankvaluechainerror() = ErrorException("all chains must have an equal length")

function _samplerenewaldidinfectionsinterventionargumenterror()
    m = "at least one of the keyword arguments `interventions` and `data` must be supplied"
    return ArgumentError(m)
end

function _samplerenewaldidinfectionsNsargumenterror()
    m = "at least one of the keyword arguments `Ns` and `data` must be supplied"
    return ArgumentError(m)
end

function _samplerenewaldidinfectionsseedmatrixargumenterror()
    m = "at least one of the keyword arguments `seedmatrix`, `data` and `observedcases` \
        must be supplied"
    return ArgumentError(m)
end

_seedwidthdimensionerror() = DimensionMismatch("width of `seedmatrix` must equal `ngroups`")
