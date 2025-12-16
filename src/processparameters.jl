# run simulation with fitted parameters 

# structs and types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
    RenewalDiDModel 

Type of the model generated from the function `renewaldid`.
"""
const RenewalDiDModel = Model{typeof(_renewaldid)}

_defaults(m::RenewalDiDModel) = m.defaults
_delaydistn(m::RenewalDiDModel) = m.args.delaydistn
_expectedseedcases(m::RenewalDiDModel) = m.args.expectedseedcases 
_generationtimefunction(m::RenewalDiDModel) = m.args.g
_interventions(m::RenewalDiDModel) = m.args.interventions 
_ngroups(m::RenewalDiDModel) = _ngroups(_interventions(m))
_ninterventions(m::RenewalDiDModel) = _ninterventions(_interventions(m))
_ns(m::RenewalDiDModel) = m.args.Ns 
_nseeds(m::RenewalDiDModel) = m.args.n_seeds
_ntimes(m::RenewalDiDModel) = _ntimes(_interventions(m))


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

julia> df = DataFrame(
       :chain => repeat([1, 2, 3]; inner=5),
       :iteration => repeat(1:5; outer=3),
       :a => rand(rng, 15)
       );

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

## quantiles of sampled outputs

"""
    quantilerenewaldidinfections(A, q=[0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975])

Calculate quantiles from an array of simulated outputs

# Arguments 
- `A`: array of outputs or a `SampledOutput`
- `q`: vector of quantiles

It is expected that an odd number of quantiles will be supplied, symmetrical around the 
    median (`0.5`)

# Examples
```jldoctest; filter=r"RenewalDiD*.*6"
julia> using StableRNGs

julia> rng = StableRNG(1000);

julia> data = RenewalDiD.testsimulation(rng);

julia> model = renewaldid(data, g_seir, RenewalDiDPriors(); mu=0.2, kappa=0.5);

julia> df = RenewalDiD.testdataframe(
       rng;
       nchains=2, niterations=5, ngroups=3, ntimes=10, nseeds=7
       );

julia> A = samplerenewaldidinfections(rng, model, df);

julia> quantilerenewaldidinfections(A.output, [0.05, 0.5, 0.95])
11×3×3 Array{Float64, 3}:
[:, :, 1] =
 0.0        0.0       0.0
 0.0        0.0       0.0
 0.0        0.0       0.0
 0.0        0.0       0.0
 0.0        0.0       0.0
 0.0        0.0       0.0
 0.131861   0.157388  0.0326495
 0.0        0.123517  0.0
 0.121198   0.759233  0.0
 0.600175   1.14889   0.0721785
 0.0508668  1.066     0.0

[:, :, 2] =
 0.0434214  0.00689775  0.0
 0.0        0.0146244   0.360681
 0.0213385  0.0413518   0.556312
 0.259588   0.0         0.312482
 0.277872   0.3432      0.146953
 0.0583196  0.410323    0.242909
 0.812611   1.02003     0.353115
 1.01233    1.9329      0.324287
 0.820946   2.58175     0.944555
 1.52041    2.68924     1.29307
 2.10001    6.06172     1.30006

[:, :, 3] =
  0.550498   0.105273  0.210392
  0.687757   0.849505  1.02244
  0.345099   1.38626   1.80516
  0.929155   1.18127   1.69581
  1.13479    0.865813  1.18232
  1.84541    3.02148   1.34952
  1.93296    2.24424   1.61729
  3.74074    4.46247   2.0209
  6.36694   12.2786    2.7994
  9.85947   13.1692    5.19129
 12.1114    25.3388    7.85569

julia> quantilerenewaldidinfections(A.output, [0.05, 0.5, 0.90]);
┌ Warning: (0.05, 0.9): other functions expect that credible intervals are symmetrical
└ @ RenewalDiD yourpath.jl:627

julia> quantilerenewaldidinfections(A.output, [0.05, 0.90]);
┌ Warning: [0.05, 0.9]: other functions expect an odd number of quantiles
└ @ RenewalDiD yourpath.jl:637
```
"""
function quantilerenewaldidinfections(
    predmodel, chains::Chains, q=[0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]; 
    #mutewarnings=nothing
)
    outputarray = Array(chains)
    ngroups = _ngroups(predmodel)
    ntimes = _ntimes(predmodel)
    ngroups * (ntimes + 1) == size(chains, 2) || throw(ErrorException("to do, $ngroups * $(ntimes + 1) == $(size(chains, 2))"))
    return _quantilerenewaldidinfections(outputarray, ngroups, ntimes, q; )
end

function quantilerenewaldidinfections(
    ::Any, array::Array{<:Any, 3}, q=[0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]; 
    #mutewarnings=nothing
)
    ngroups = size(array, 3)
    ntimes = size(array, 2) - 1
    return _quantilerenewaldidinfections(array, ngroups, ntimes, q; )
end

function _quantilerenewaldidinfections(outputarray::Array{<:Any, 2}, ngroups, ntimes, q; )
    outputquantiles = zeros(ntimes + 1, ngroups, length(q))

    i = 1
    for g in 1:ngroups, t in 1:(ntimes + 1)
        outputquantiles[t, g, :] .= max.(0, quantile(skipmissing(outputarray[:, i]), q))
        i += 1
    end

    return outputquantiles
end

function _quantilerenewaldidinfections(outputarray::Array{<:Any, 3}, ngroups, ntimes, q; )
    outputquantiles = zeros(ntimes + 1, ngroups, length(q))
    for g in 1:ngroups, t in 1:(ntimes + 1)
        outputquantiles[t, g, :] .= max.(0, quantile(skipmissing(outputarray[:, t, g]), q))
    end
    return outputquantiles
end

## list the intervention times

# other versions of this function are in `interventionmatrix.jl`
_interventionstarttimes(M::AbstractMatrix, i) = findfirst(x -> x == 1, M[:, i])

## versions of gamma and theta vector functions taking a `DataFrame` of fitted parameters

function _gammavec(df::DataFrame, i, ngroups)
    return _gammavec(  # call version in `fittingparameters.jl`
        [getproperty(df, Symbol("gammas_raw[$j]"))[i] for j in 1:(ngroups - 1)], 
        df.sigma_gamma[i]
    )
end 

function _thetavec(df::DataFrame, i, ntimes; thetainterval=automatic)
    nthetas = _nthetas(ntimes, thetainterval)
    return _thetavec(  # call version in `fittingparameters.jl`
        [getproperty(df, Symbol("thetas_raw[$j]"))[i] for j in 1:nthetas], 
        df.sigma_theta[i],
        ntimes;
        thetainterval,
    )
end  

function _tauvec(df, i, ninterventions)
    return [getproperty(df, Symbol("logtau[$k]"))[i] for k in 1:ninterventions]
end  

function _mxmatrix(df, i, ngroups, ntimes, n_seeds)
    totaltimes = ntimes + n_seeds
    return [getproperty(df, Symbol("M_x[$t, $j]"))[i] for t in 1:totaltimes, j in 1:ngroups]
end

predictedR_0(model::RenewalDiDModel, chain::Chains) = predictedR_0(model, DataFrame(chain))
predictedR_0(model::RenewalDiDModel, chaindf::DataFrame) = _modelpredictedR_0(model, chaindf)

function _modelpredictedR_0(model, chaindf)
    nsamples = size(chaindf, 1)
    ntimes = _ntimes(model)
    ngroups = _ngroups(model)
    ninterventions = _ninterventions(model)
    return _modelpredictedR_0(model, chaindf, nsamples, ntimes, ngroups, ninterventions)
end

function _modelpredictedR_0(model, chaindf, nsamples, ntimes, ngroups, ninterventions)
    R0s = zeros(nsamples, ntimes, ngroups)
    _modelpredictedR_0!(R0s, model, chaindf, nsamples, ntimes, ngroups, ninterventions)
    return R0s
end

function _modelpredictedR_0!(R0s, model, chaindf, nsamples, ntimes, ngroups, ninterventions)
    for i in 1:nsamples 
        R0s[i, :, :] .= _modelpredictedR_0iteration(
            model, chaindf, ntimes, ngroups, ninterventions, i
        )
    end
    return nothing
end

function _modelpredictedR_0iteration(model, chaindf, ntimes, ngroups, ninterventions, i)
    alpha = chaindf.alpha[i]
    gammavec = _gammavec(chaindf, i, ngroups)
    thetavec = _thetavec(chaindf, i, ntimes; thetainterval=model.args.thetainterval)
    tau = _tauvec(chaindf, i, ninterventions)
    return exp.(
        _predictedlogR_0(alpha, gammavec, thetavec, tau, model.args.interventions)
    )
end

function predictedcases(predmodel, chain; kwargs...)
    return predictedcases(default_rng(), predmodel, DataFrame(chain); kwargs...)
end

function predictedcases(rng::AbstractRNG, predmodel, chain::Chains; kwargs...)
    return predictedcases(rng, predmodel, DataFrame(chain); kwargs...)
end

function predictedcases(rng::AbstractRNG, predmodel, chaindf::DataFrame; kwargs...)
    nsamples = size(chaindf, 1)
    ntimes = _ntimes(predmodel)
    ngroups = _ngroups(predmodel)
    ninterventions = _ninterventions(predmodel)
    nseeds = _nseeds(predmodel)
    Ns = _ns(predmodel)
    return _predictedcases(
        rng, predmodel, chaindf, nsamples, ntimes, ngroups, ninterventions, Ns, nseeds; 
        kwargs...
    )
end

function _predictedcases(
    rng::AbstractRNG, 
    predmodel, 
    chaindf, 
    nsamples, 
    ntimes, 
    ngroups, 
    ninterventions, 
    Ns, 
    nseeds; 
    kwargs...
)
    outputcases = Array{Union{Float64, Missing}}(undef, nsamples, ntimes + 1, ngroups)
    R0s = _modelpredictedR_0(predmodel, chaindf, nsamples, ntimes, ngroups, ninterventions)
    _predictedcases!(
        outputcases, rng, predmodel, R0s, chaindf, nsamples, ntimes, ngroups, Ns, nseeds; 
        kwargs...
    ) 
    return outputcases
end

function _predictedcases!(
    outputcases, 
    rng::AbstractRNG, 
    predmodel, 
    R0s, 
    chaindf, 
    nsamples, 
    ntimes, 
    ngroups,
    Ns, 
    nseeds; 
    kwargs...
)
    predictedinfections = Array{Float64}(undef, ntimes + nseeds, ngroups)

    for i in 1:nsamples 
        _predictedcasesiteration!(
            outputcases, 
            predictedinfections, 
            rng, 
            predmodel, 
            R0s, 
            chaindf, 
            ntimes, 
            ngroups, 
            Ns, 
            nseeds, 
            i; 
            kwargs...
        )
    end
    return nothing
end

function _predictedcasesiteration!(
    outputcases, 
    predictedinfections, 
    rng::AbstractRNG, 
    predmodel, 
    R0s, 
    chaindf, 
    ntimes, 
    ngroups, 
    Ns, 
    nseeds, 
    i; 
    maxdelay=automatic, smoothmaxepsilon=0.1, kwargs...
)
    psi = getproperty(chaindf, :psi)[i]
    negbinom_r = getproperty(chaindf, :negbinom_r)[i]

    for j in 1:ngroups 
        s = one(Float64)  # track `s` here rather than complex numbers 
        for t in 1:(nseeds + ntimes) 
            if t <= nseeds 
                predictedinfections[t, j] = predmodel.args.expectedseedcases[t, j] * s
            else
                predictedinfections[t, j] = 
                    R0s[i, (t - nseeds), j] * 
                    sum(
                        [
                            predictedinfections[x, j] * predmodel.args.g(t - x; kwargs...) 
                            for x in 1:(t - 1)
                        ]
                    ) * 
                    s
            end

            s = RenewalDiD._track_s(s, predictedinfections, Ns, t, j; epsilon=smoothmaxepsilon)
        end
    end

    # delay between infection and detection
    delayedinfections = _delayedinfections(
        Float64, predictedinfections, predmodel.args.delaydistn, ngroups, ntimes, nseeds, maxdelay
    )

    # negative binomail likelihood 
    negbinom_p = _negbinom_p(
        negbinom_r, delayedinfections[nseeds:(nseeds + ntimes), :] .* psi
    )

    if isnan(max(negbinom_r, maximum(negbinom_p))) 
        outputcases[i, :, :] .= [missing for _ in 0:ntimes, _ in 1:ngroups]
    else
        outputcases[i, :, :] .= rand.(rng, NegativeBinomial.(negbinom_r, negbinom_p))
    end

    return nothing
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
