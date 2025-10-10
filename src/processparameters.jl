# run simulation with fitted parameters 

# structs and types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
    SampledOutput{T, N} 

Contains the number of reported infections and the basic reproduction ratio output from a 
    model.

# Fields
- `output::Array{T, N}`: numbers of infections
- `R0s::Array{T, N}`: time-varying basic reproduction ratio
"""
struct SampledOutput{T, N} 
    output::Array{T, N}
    R0s::Array{T, N}

    function SampledOutput(output::Array{T, N}, R0s::Array{T, N}) where {T, N}
        return new{T, N}(output, R0s)
    end
end

Base.hash(x::SampledOutput, h::UInt64) = hash(x.R0s, hash(x.output, hash(:SampledOutput, h)))

function Base.:(==)(a::SampledOutput, b::SampledOutput)
    c1 = a.output == b.output 
    c2 = a.R0s == b.R0s 
    if c1 === false || c2 === false 
        return false 
    elseif ismissing(c1) || ismissing(c2)
        return missing 
    end
    return true 
end

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
        df.sigma_theta[i],
        ntimes,
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
    samplerenewaldidinfections([rng], model, fittedparameters[, indexes; <keyword arguments>])

Generate expected outcomes from model using samples from distributions.    

# Arguments 
- `rng::AbstractRNG=default_rng()`: random number generator
- `model::RenewalDiDModel`: model being used 
- `fittedparameters`: can be a `DataFrame`, `Chains` or `Tuple{DataFrame, <:Chains}`
- `indexes::Union{<:AbstractVector{<:Integer}, <:Integer}`: samples to use, defaults to all 
    samples

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(1000);

julia> data = RenewalDiD.testsimulation(rng);

julia> model = renewaldid(data, g_seir, RenewalDiDPriors(); mu=0.2, kappa=0.5);

julia> df = RenewalDiD.testdataframe(
       rng;
       nchains=2, niterations=5, ngroups=3, ntimes=10, nseeds=7
       );

julia> samplerenewaldidinfections(rng, model, df, 1)
SampledOutput{Float64, 2}([0.07491152793594207 0.0 0.0; 0.0 0.1577299327763031 \
    0.2629712097378153; … ; 1.2264839609837064 2.3251829893883746 0.16039663151955016; \
    0.740412457261556 3.8074225038310803 0.0], [1.4090694615560517 1.4090694615560517 \
    0.0656949120738175; 1.446863906052742 1.446863906052742 0.10348935657050783; … ; \
    1.7114250175295744 2.283489655473959 0.9401151059917245; 1.7492194620262647 \
    2.321284099970649 0.9779095504884148])
```
"""
function samplerenewaldidinfections(
    model::RenewalDiDModel, fittedparameters, indexes=automatic; 
    kwargs...
)
    return samplerenewaldidinfections(
        default_rng(), model, fittedparameters, indexes; 
        kwargs...
    )
end

function samplerenewaldidinfections(
    rng::AbstractRNG, 
    model::RenewalDiDModel, 
    fittedparameters::Tuple{DataFrame, <:Chains}, 
    indexes=automatic; 
    kwargs...
)
    return samplerenewaldidinfections(rng, model, fittedparameters[1], indexes; kwargs...)
end

function samplerenewaldidinfections(
    rng::AbstractRNG, model::RenewalDiDModel, fittedparameters::Chains, indexes=automatic; 
    kwargs...
)
    return samplerenewaldidinfections(
        rng, model, DataFrame(fittedparameters), indexes; 
        kwargs...
    )
end

function samplerenewaldidinfections(
    rng::AbstractRNG, model::RenewalDiDModel, fittedparameters::DataFrame, indexes=automatic; 
    repeatsamples=nothing, kwargs...
)
    return _samplerenewaldidinfections(
        rng, model, fittedparameters, indexes, repeatsamples; 
        kwargs...
    )
end

function samplerenewaldidinfections(
    ::AbstractRNG, ::RenewalDiDModel, fittedparameters, x=nothing; 
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
    ngroups = _ngroups(model)
    ntimes = _ntimes(model)
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
    ngroups = _ngroups(model)
    ntimes = _ntimes(model)
    _samplerenewaldidinfectionsassertions(
        fittedparameters, _expectedseedcases(model), i, ngroups, ntimes
    )
    alpha = fittedparameters.alpha[i]
    gammavec = _gammavec(fittedparameters, i, ngroups)
    thetavec = _thetavec(fittedparameters, i, ntimes)
    tau = _tauvec(fittedparameters, i, _ninterventions(model))
    psi = fittedparameters.psi[i]
    minsigma2 = fittedparameters.minsigma2[i]
    M_x = _mxmatrix(fittedparameters, i, ngroups, ntimes, _nseeds(model))
    R0s .= _predictedlogR_0(alpha, gammavec, thetavec, tau, _interventions(model))
    T = Complex{typeof(R0s[1, 1])}
    predictedinfections = _infectionsmatrix(T, R0s, _nseeds(model))
    _infections!(
        _generationtimefunction(model), 
        predictedinfections, 
        M_x, 
        R0s, 
        _expectedseedcases(model), 
        _ns(model), 
        _nseeds(model); 
        _defaults(model)...
    )
    delayedinfections = _delayedinfections(
        T, predictedinfections, _delaydistn(model), ngroups, ntimes, _nseeds(model)
    )
    np = real.(delayedinfections[_nseeds(model):(_nseeds(model) + ntimes), :]) .* psi
    isnan(maximum(np)) && return _halt_samplerenewaldidinfections(i)  # exit early 
    output .= _predictobservedinfections(rng, np, psi, minsigma2)
    return nothing
end

function _predictobservedinfections(rng, np::AbstractMatrix, psi, minsigma2)
    return [
        __predictobservedinfections(rng, np[t, g], psi, minsigma2) 
        for t in axes(np, 1), g in axes(np, 2)
    ]
end

function __predictobservedinfections(rng, np::Number, psi, minsigma2)
    v = rand(rng, Normal(np, sqrt(np * (1 - psi) + minsigma2)))
    v > 0 || return zero(v)
    return v
end

function _halt_samplerenewaldidinfections(i) 
    @warn "Parameters in row $i give `NaN` predicted infections"
    return nothing 
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
    SO::SampledOutput, q=[0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]; 
    mutewarnings=nothing
)
    return SampledOutput(
        quantilerenewaldidinfections(SO.output, q; mutewarnings),
        quantilerenewaldidinfections(SO.R0s, q; mutewarnings)
    )
end

function quantilerenewaldidinfections(
    A::AbstractArray, q=[0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]; 
    mutewarnings=nothing
)
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

#function map_DataFrame(result::ModeResult)
function map_DataFrame(result)
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
