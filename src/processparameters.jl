# run simulation with fitted parameters 

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

function _mxmatrix(df, i, ngroups, ntimes, n_seeds)
    totaltimes = ntimes + n_seeds
    return [getproperty(df, Symbol("M_x[$t, $j]"))[i] for t in 1:totaltimes, j in 1:ngroups]
end

## take samples from DataFrame of fitted parameters and generate expected outcomes

function samplerenewaldidinfections(
    g, df::DataFrame, indexes::AbstractVector{<:Integer}=axes(df, 1); 
    kwargs...
)
    interventions, Ns, seedmatrix, ngroups, ntimes, n_seeds, kws = __sra(; kwargs...)
    infn = zeros(Float64, ntimes + n_seeds, ngroups, length(indexes))
    _samplerenewaldidinfectionsassertions(df, seedmatrix, indexes, ngroups, ntimes)
    for (r, j) in enumerate(indexes)
        _samplerenewaldidinfections!(
            g, 
            (@view infn[:, :, r]), 
            df, 
            interventions, 
            Ns, 
            seedmatrix, 
            j, 
            ngroups, 
            ntimes, 
            n_seeds;
            kws...
        )
    end
    return infn[n_seeds:n_seeds+ntimes, :, :]
end

function samplerenewaldidinfections(g, df::DataFrame, i::Integer; kwargs...)
    interventions, Ns, seedmatrix, ngroups, ntimes, n_seeds, kws = __sra(; kwargs...)
    infn = zeros(Float64, ntimes + n_seeds, ngroups)
    _samplerenewaldidinfectionsassertions(df, seedmatrix, i, ngroups, ntimes)
    _samplerenewaldidinfections!(
        g, infn, df, interventions, Ns, seedmatrix, i, ngroups, ntimes, n_seeds;
        kws...
    )
    return infn[n_seeds:n_seeds+ntimes, :]
end

function samplerenewaldidinfections!(g, infn, df::DataFrame, i; kwargs...)
    interventions, Ns, seedmatrix, ngroups, ntimes, n_seeds, kws = __sra(; kwargs...)
    _samplerenewaldidinfectionsassertions(df, seedmatrix, i, ngroups, ntimes)
    _samplerenewaldidinfections!(
        g, infn, df, interventions, Ns, seedmatrix, i, ngroups, ntimes, n_seeds;
        kws...
    )
    return nothing
end

function _samplerenewaldidinfections!(
    g::Vector, infn, df, interventions, Ns, seedmatrix, i, ngroups, ntimes, n_seeds;
    kwargs...
)
    return _samplerenewaldidinfections!(
        generationtime, infn, df, interventions, Ns, seedmatrix, i, ngroups, ntimes, n_seeds; 
        vec=g, kwargs...
    )
end

function _samplerenewaldidinfections!(
    g::_Useablegenerationfunctions, 
    infn, 
    df, 
    interventions, 
    Ns, 
    seedmatrix, 
    i, 
    ngroups, 
    ntimes, 
    n_seeds; 
    kwargs...
)
    alpha = df.alpha[i]
    gammavec = _gammavec(df, i, ngroups)
    thetavec = _thetavec(df, i, ntimes)
    tau = df.tau[i]
    logR_0 = _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)
    M_x = _mxmatrix(df, i, ngroups, ntimes, n_seeds)
    _infections!(
        g, 
        infn,
        M_x, 
        logR_0, 
        seedmatrix, 
        Ns, 
        n_seeds; 
        kwargs...
    )
    return nothing
end

## quantiles of sampled outputs

function quantilerenewaldidinfections(A, q; mutewarnings=nothing)
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


# select parameters from keyword arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

__sra(; kwargs...) = _samplerenewaldidinfectionarguments(; kwargs...)

function _samplerenewaldidinfectionarguments( ;
    data=nothing,
    doubletime=automatic,
    interventions=nothing, 
    n_seeds=nothing,
    ngroups=automatic, 
    ntimes=automatic,
    Ns=nothing, 
    observedcases=nothing,
    sampletime=automatic,
    seedcasesminvalue=0.5,
    seedmatrix=automatic, 
    kwargs...
)
    return __samplerenewaldidinfectionarguments( 
        data,
        doubletime,
        interventions, 
        n_seeds,
        ngroups, 
        ntimes,
        Ns, 
        observedcases,
        sampletime,
        seedcasesminvalue,
        seedmatrix;
        kwargs...
    )
end

function __samplerenewaldidinfectionarguments( 
    inputdata,
    doubletime,
    inputinterventions, 
    inputn_seeds,
    inputngroups, 
    inputntimes,
    inputNs, 
    inputobservedcases,
    sampletime,
    seedcasesminvalue,
    inputseedmatrix;
    kwargs...
)
    interventions = _samplerenewaldidinfectionsinterventions(inputdata, inputinterventions)
    Ns = _samplerenewaldidinfectionsNs(inputdata, inputNs)
    seedmatrix = _samplerenewaldidinfectionsseedmatrix(
        inputdata, inputseedmatrix, inputobservedcases, inputn_seeds; 
        minvalue=seedcasesminvalue, doubletime, sampletime
    )
    ngroups = _samplerenewaldidinfectionsngroups(inputngroups, interventions)
    ntimes = _samplerenewaldidinfectionsntimes(inputntimes, interventions)
    n_seeds = _samplerenewaldidinfectionsnseeds(inputn_seeds, seedmatrix)
    kws = kwargs
    return (interventions, Ns, seedmatrix, ngroups, ntimes, n_seeds, kws)
end

function _samplerenewaldidinfectionsinterventions(::Any, inputinterventions::AbstractMatrix)
    return inputinterventions
end

function _samplerenewaldidinfectionsinterventions(inputdata::Dict, ::Nothing)
    return inputdata[:interventions]
end

function _samplerenewaldidinfectionsinterventions(::Nothing, ::Nothing)
    throw(_samplerenewaldidinfectionsinterventionargumenterror())
    return nothing
end

_samplerenewaldidinfectionsNs(::Any, Ns::AbstractVector) = Ns
_samplerenewaldidinfectionsNs(inputdata::Dict, ::Nothing) = inputdata[:Ns]

function _samplerenewaldidinfectionsNs(::Nothing, ::Nothing)
    throw(_samplerenewaldidinfectionsNsargumenterror())
    return nothing
end

function _samplerenewaldidinfectionsseedmatrix(
    ::Any, seedmatrix::AbstractMatrix, ::Any, ::Any; 
    kwargs...
)
    return seedmatrix
end

function _samplerenewaldidinfectionsseedmatrix( 
    ::Any, ::Automatic, inputobservedcases::AbstractMatrix, n_seeds::Integer; 
    kwargs...
)
    return _expectedseedcases(inputobservedcases, n_seeds; kwargs...)
end 

function _samplerenewaldidinfectionsseedmatrix( 
    inputdata::Dict, ::Automatic, ::Any, n_seeds::Any; 
    kwargs...
)
    return _samplerenewaldidinfectionsseedmatrix(
        nothing, automatic, inputdata[:observedcases], n_seeds; 
        kwargs...
    )
end

function _samplerenewaldidinfectionsseedmatrix(
    ::Nothing, ::Automatic, ::Nothing, ::Any; 
    kwargs...
)
    throw(_samplerenewaldidinfectionsseedmatrixargumenterror())
    return nothing
end

function _samplerenewaldidinfectionsseedmatrix( 
    inputdata::Any, ::Automatic, inputobservedcases::AbstractMatrix, ::Nothing; 
    kwargs...
)
    # if n_seeds is not supplied, use the default from `renewaldid`
    return _samplerenewaldidinfectionsseedmatrix(
        inputdata, automatic, inputobservedcases, 7; 
        kwargs...
    )
end

_samplerenewaldidinfectionsngroups(inputngroups::Integer, ::Any) = inputngroups
_samplerenewaldidinfectionsngroups(::Automatic, interventions) = _ngroups(interventions)
_samplerenewaldidinfectionsntimes(inputntimes::Integer, ::Any) = inputntimes
_samplerenewaldidinfectionsntimes(::Automatic, interventions) = _ntimes(interventions)
_samplerenewaldidinfectionsnseeds(inputn_seeds::Integer, ::Any) = inputn_seeds
_samplerenewaldidinfectionsnseeds(::Nothing, seedmatrix) = _ntimes(seedmatrix)


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

function _quantilerenewaldidinfectionswarningset(A::AbstractArray{T, 3}, q, ::Nothing) where T
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
