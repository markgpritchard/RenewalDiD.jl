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

function _observedsigmamatrix(df, i, ngroups, ntimes)
    return return [getproperty(df, Symbol("predictobservedinfectionssigmamatrix[$t, $j]"))[i] for t in 1:(ntimes + 1), j in 1:ngroups]
end

## take samples from DataFrame of fitted parameters and generate expected outcomes

function samplerenewaldidinfections(
    g, df::DataFrame, data::AbstractRenewalDiDData, indexes::AbstractVector{<:Integer}=axes(df, 1);
    kwargs...
)
    ngroups = _ngroups(data.interventions)
    ntimes = _ntimes(data.interventions)
    n_seeds = size(data.exptdseedcases, 1)
    output = zeros(ntimes + 1, ngroups, length(indexes))
    for (r, j) in enumerate(indexes)
        _samplerenewaldidinfections!(
            g, (@view output[:, :, r]), df, data, j, ngroups, ntimes, n_seeds; 
            kwargs...
        )
    end
    return output 
end

function samplerenewaldidinfections(
    g, df::DataFrame, data::AbstractRenewalDiDData, i::Integer;
    kwargs...
)
    ngroups = _ngroups(data.interventions)
    ntimes = _ntimes(data.interventions)
    n_seeds = size(data.exptdseedcases, 1)
    output = zeros(ntimes + 1, ngroups)
    _samplerenewaldidinfections!(
            g, output, df, data, i, ngroups, ntimes, n_seeds; 
            kwargs...
        )
    return output 
end

function _samplerenewaldidinfections!(
    g, output::AbstractArray, df, data::AbstractRenewalDiDData, i, ngroups, ntimes, n_seeds; 
    kwargs...
)
    return _samplerenewaldidinfections!(
        generationtime, output, df, data, i, ngroups, ntimes, n_seeds; 
        func=g,  # `generationtime` accepts function or vector from the keyword `func`
        kwargs...
    )
end

function _samplerenewaldidinfections!(
    g::_Useablegenerationfunctions, 
    output::AbstractArray, 
    df, 
    data::RenewalDiDData, 
    i, 
    ngroups, 
    ntimes, 
    n_seeds; 
    kwargs...
)
    Ns = data.Ns  
    return _samplerenewaldidinfections!(
        g, output::AbstractArray, df, data, i, ngroups, ntimes, n_seeds, Ns; 
        kwargs...
    )
end

function _samplerenewaldidinfections!(
    g::_Useablegenerationfunctions, 
    output::AbstractArray, 
    df, 
    data::RenewalDiDDataUnlimitedPopn, 
    i, 
    ngroups, 
    ntimes, 
    n_seeds; 
    kwargs...
)
    return _samplerenewaldidinfections!(
        g, output::AbstractArray, df, data, i, ngroups, ntimes, n_seeds, nothing; 
        kwargs...
    )
end

function _samplerenewaldidinfections!(
    g::_Useablegenerationfunctions, 
    output::AbstractArray, 
    df, 
    data::AbstractRenewalDiDData, 
    i, 
    ngroups, 
    ntimes, 
    n_seeds, 
    Ns; 
    kwargs...
)
    alpha = df.alpha[i]
    gammavec = _gammavec(df, i, ngroups)
    thetavec = _thetavec(df, i, ntimes)
    tau = df.tau[i]
    psi = df.psi[i]
    M_x = _mxmatrix(df, i, ngroups, ntimes, n_seeds)
    predictobservedinfectionssigmamatrix = _observedsigmamatrix(df, i, ngroups, ntimes)

    predictedlogR_0 = _predictedlogR_0(alpha, gammavec, thetavec, tau, data.interventions)
    T = Complex{typeof(predictedlogR_0[1, 1])}
    predictedinfections = _infectionsmatrix(T, predictedlogR_0, n_seeds)
    _infections!(
        g, predictedinfections, M_x, predictedlogR_0, data.exptdseedcases, Ns, n_seeds; 
        kwargs...
    )
    np = real.(predictedinfections[n_seeds:n_seeds+ntimes, :]) .* psi
    isnan(maximum(np)) && return _halt_samplerenewaldidinfections(i)  # exit early 
    predictobservedinfections = max.(
        0, 
        np .+ predictobservedinfectionssigmamatrix .* np .* (1 - psi)
    )
    output .= predictobservedinfections
    return nothing
end

function _halt_samplerenewaldidinfections(i) 
    @warn "Parameters in row $i give `NaN` predicted infections"
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
