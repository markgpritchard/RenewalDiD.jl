# run simulation with fitted parameters 

## Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### number unique

nunique(vec) = length(unique(vec))

### data for tracerankplot 

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



### versions of gamma and theta vector functions 

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

### take samples from DataFrame of fitted parameters

function samplerenewaldidinfections(
    g, df, i; 
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
    return _samplerenewaldidinfections(
        g, 
        df, 
        i, 
        data, 
        interventions, 
        ngroups, 
        ntimes, 
        Ns, 
        observedcases, 
        seedmatrix, 
        n_seeds; 
        doubletime, sampletime, seedcasesminvalue,
        kwargs...
    )
end

function _samplerenewaldidinfections(
    g, 
    df, 
    i, 
    inputdata, 
    inputinterventions, 
    inputngroups, 
    inputntimes, 
    inputNs, 
    inputobservedcases, 
    inputseedmatrix, 
    n_seeds; 
    doubletime, sampletime, seedcasesminvalue,
    kwargs...
)
    interventions = _samplerenewaldidinfectionsinterventions(inputdata, inputinterventions)
    Ns = _samplerenewaldidinfectionsNs(inputdata, inputNs)
    seedmatrix = _samplerenewaldidinfectionsseedmatrix(
        inputdata, inputseedmatrix, inputobservedcases, n_seeds; 
        minvalue=seedcasesminvalue, doubletime, sampletime
    )
    ngroups = _samplerenewaldidinfectionsngroups(inputngroups, interventions)
    ntimes = _samplerenewaldidinfectionsntimes(inputntimes, interventions)
    _samplerenewaldidinfectionsassertions(df, seedmatrix, i, ngroups, ntimes)
    return _samplerenewaldidinfections(
        g, df, interventions, Ns, seedmatrix, i, ngroups, ntimes;
        kwargs...
    )
end

function _samplerenewaldidinfections(
    g::Vector, df, interventions, Ns, seedmatrix, i, ngroups, ntimes;
    kwargs...
)
    return _samplerenewaldidinfections(
        generationtime, df, interventions, Ns, seedmatrix, i, ngroups, ntimes; 
        vec=g, kwargs...
        )
end

function _samplerenewaldidinfections(
    g::_Useablegenerationfunctions, df, interventions, Ns, seedmatrix, i, ngroups, ntimes; 
    kwargs...
)
    n_seeds = _ntimes(seedmatrix)
    alpha = df.alpha[i]
    gammavec = _gammavec(df, i, ngroups)
    thetavec = _thetavec(df, i, ntimes)
    tau = df.tau[i]
    logR_0 = _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)
    infn = _infections(
        Float64, 
        g, 
        zeros(ntimes + n_seeds, ngroups), 
        logR_0, 
        seedmatrix, 
        Ns, 
        n_seeds; 
        kwargs...
    )
    return infn[n_seeds:n_seeds+ntimes, :]
end

#### select parameters from keyword arguments

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
    _expectedseedcases(observedcases, n_seeds; kwargs...)
end 

function _samplerenewaldidinfectionsseedmatrix( 
    inputdata::Dict, ::Automatic, ::Any, n_seeds::Any; 
    kwargs...
)
    _samplerenewaldidinfectionsseedmatrix(
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
    ::Any, ::Automatic, ::AbstractMatrix, ::Nothing; 
    kwargs...
)
    throw(_samplerenewaldidinfectionsseedmatrixnseedserror())
    return nothing
end

_samplerenewaldidinfectionsngroups(inputngroups::Integer, ::Any) = inputngroups
_samplerenewaldidinfectionsngroups(::Automatic, interventions) = _ngroups(interventions)
_samplerenewaldidinfectionsntimes(inputntimes::Integer, ::Any) = inputntimes
_samplerenewaldidinfectionsntimes(::Automatic, interventions) = _ntimes(interventions)


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
    i <= size(df, 1) || throw(BoundsError(df, [i, 1]))
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

function _samplerenewaldidinfectionsseedmatrixnseedserror()
    m = "if the keyword argument `seedmatrix` is not supplied then `n_seeds` must be"
    return ArgumentError(m)
end

_seedwidthdimensionerror() = DimensionMismatch("width of `seedmatrix` must equal `ngroups`")
