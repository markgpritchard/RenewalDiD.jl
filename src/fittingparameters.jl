# Use Turing.jl to fit parameters to a dataset

# Structs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## data passed to `renewaldid`

abstract type AbstractRenewalDiDData end

@auto_hash_equals struct RenewalDiDData{S, T} <: AbstractRenewalDiDData
    observedcases::Matrix{S}
    interventions::T
    Ns::Vector{Int}
    exptdseedcases::Matrix{Float64}

    function RenewalDiDData(
        observedcases::Matrix{S}, 
        interventions::T, 
        Ns::Vector{<:Number}, 
        exptdseedcases::Matrix{<:Number}
    ) where {S <:Number, T <:AbstractArray}
        _diddataassertions(observedcases, interventions, Ns, exptdseedcases)
        return new{S, T}(
            observedcases, 
            interventions, 
            convert.(Int, Ns), 
            convert.(Float64, exptdseedcases)
        )
    end
end

function RenewalDiDData( ; 
    observedcases, 
    interventions, 
    Ns, 
    exptdseedcases=nothing,
    n_seeds=DEFAULT_SEEDMATRIX_HEIGHT, 
    doubletime=automatic, 
    sampletime=automatic, 
    minvalue=DEFAULT_SEEDMATRIX_MINVALUE,
)
    newexptdseedcases = _expectedseedcasesifneeded(
        exptdseedcases, observedcases, n_seeds; 
        doubletime, sampletime, minvalue
    )
    return RenewalDiDData(observedcases, interventions, Ns, newexptdseedcases)
end

@auto_hash_equals struct RenewalDiDDataUnlimitedPopn{S, T} <: AbstractRenewalDiDData
    observedcases::Matrix{S}
    interventions::T
    exptdseedcases::Matrix{Float64}

    function RenewalDiDDataUnlimitedPopn(
        observedcases::Matrix{S}, 
        interventions::T, 
        exptdseedcases::Matrix{<:Number}
    ) where {S <:Number, T <:AbstractArray}
        _diddataassertions(observedcases, interventions, exptdseedcases)
        return new{S, T}(
            observedcases, 
            interventions, 
            convert.(Float64, exptdseedcases)
        )
    end
end

function RenewalDiDDataUnlimitedPopn( ; 
    observedcases, 
    interventions, 
    Ns=nothing,
    exptdseedcases=nothing,
    n_seeds=DEFAULT_SEEDMATRIX_HEIGHT, 
    doubletime=automatic, 
    sampletime=automatic, 
    minvalue=DEFAULT_SEEDMATRIX_MINVALUE,
)
    return _RenewalDiDDataUnlimitedPopn(
        observedcases, 
        interventions, 
        Ns,
        exptdseedcases,
        n_seeds, 
        doubletime, 
        sampletime, 
        minvalue,
    )
end

function _RenewalDiDDataUnlimitedPopn(
    observedcases, 
    interventions, 
    Ns,
    exptdseedcases,
    n_seeds, 
    doubletime, 
    sampletime, 
    minvalue,
)
    @warn "$Ns: keyword argument `Ns` not used in constructing `RenewalDiDDataUnlimitedPopn`"
    return _RenewalDiDDataUnlimitedPopn(
        observedcases, 
        interventions, 
        nothing,
        exptdseedcases,
        n_seeds, 
        doubletime, 
        sampletime, 
        minvalue,
    )
end

function _RenewalDiDDataUnlimitedPopn(
    observedcases, 
    interventions, 
    ::Nothing,
    exptdseedcases,
    n_seeds, 
    doubletime, 
    sampletime, 
    minvalue,
)
    newexptdseedcases = _expectedseedcasesifneeded(
        exptdseedcases, observedcases, n_seeds; 
        doubletime, sampletime, minvalue
    )
    return RenewalDiDDataUnlimitedPopn(observedcases, interventions, newexptdseedcases)
end

function Base.getproperty(obj::RenewalDiDDataUnlimitedPopn, name::Symbol)
    if name === :Ns 
        return nothing 
    else 
        return getfield(obj, name)
    end
end

function Base.show(io::IO, ::MIME"text/plain", d::AbstractRenewalDiDData)
    _showtitle(io, d)
    print(io, "\n observedcases:  ")
    show(io, d.observedcases)
    print(io, "\n interventions:  ")
    show(io, d.interventions)
    print(io, "\n Ns:             ")
    _showns(io, d)
    print(io, "\n exptdseedcases: ")
    show(io, d.exptdseedcases)
    return nothing
end

_showtitle(io, ::RenewalDiDData{S, T}) where {S, T} = print(io, "RenewalDiDData{$S, $T}")

function _showtitle(io, ::RenewalDiDDataUnlimitedPopn{S, T}) where {S, T}
    print(io, "RenewalDiDDataUnlimitedPopn{$S, $T}")
    return nothing
end

_showns(io, d::RenewalDiDData) = show(io, d.Ns)
_showns(io, ::RenewalDiDDataUnlimitedPopn) = print(io, "unlimited")

## struct of prior distributions 

@kwdef struct RenewalDiDPriors{Q, R, S, T, U, V, W, X, Y}
    # attempting to fit `mu_delay` and `sigma_delay` gives NaN gradients so currently use
    # constants. `renewaldid` will throw a MethodError if these are distributions
    alphaprior::Q=Normal(0, 1)
    M_xprior::R=truncated(Normal(0, 1); lower=-1)
    mu_delayprior::S=log(2)  #mu_delayprior::S=Normal(0, 1)
    sigma_delayprior::T=log(5)  #sigma_delayprior::T=Exponential(1)
    sigma_gammaprior::U=Exponential(1)
    sigma_thetaprior::V=Exponential(1)
    tauprior::W=Normal(0, 1)
    psiprior::X=Beta(1, 1)
    omegaprior::Y=0
end

function Base.show(
    io::IO, ::MIME"text/plain", p::RenewalDiDPriors{Q, R, S, T, U, V, W, X, Y}
) where {Q, R, S, T, U, V, W, X, Y}
    print(
        io,
        "RenewalDiDPriors{$Q, $R, $S, $T, $U, $V, $W, $X, $Y}",
        "\n alphaprior:       ", p.alphaprior,
        "\n M_xprior:         ", p.M_xprior,
        "\n mu_delayprior:    ", p.mu_delayprior,
        "\n sigma_delayprior: ", p.sigma_delayprior,
        "\n sigma_gammaprior: ", p.sigma_gammaprior,
        "\n sigma_thetaprior: ", p.sigma_thetaprior,
        "\n tauprior:         ", p.tauprior,
        "\n psiprior:         ", p.psiprior,
        "\n omegaprior:       ", p.omegaprior,
    )
end


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

## Functions called by `_renewaldid`

logsumexp(args...; alpha) = log(sum([exp(alpha * x) for x in args])) / alpha

_ntimes(M::AbstractMatrix) = size(M, 1)
_ngroups(M::AbstractMatrix) = size(M, 2)

# the mean of all gamma values is zero, ensured by setting the final value in the vector as 
# `-sum` of other values
_gammavec(gammas_raw, sigma_gamma) = [gammas_raw; -sum(gammas_raw)] .* sigma_gamma
# a second version of this function is given in `processparameters.jl`

# The value of theta at time 0 is fixed to 0. `thetas_raw` is a vector representing how much 
# each subsequent theta differs from the previous one as a multiple of the standard 
# deviation. `sigma_theta` is the standard deviation. The output of this function is the 
# time-varying `theta` at each time as a cumulative random walk. 
function _thetavec(thetas_raw::AbstractVector{S}, sigma_theta::T) where {S, T}
    theta_0 = zero(S) * zero(T)
    return cumsum([theta_0; thetas_raw .* sigma_theta])
end  # a second version of this function is given in `processparameters.jl`

function _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)
    _predictedlogR_0assertions(gammavec, thetavec, interventions)
    logR_0 = [
        alpha + gammavec[g] + thetavec[t] + tau * interventions[t, g] 
        for t in axes(interventions, 1), g in axes(interventions, 2)
    ]
    return logR_0
end

function expectedseedcases(
    observedcases, n_seeds=DEFAULT_SEEDMATRIX_HEIGHT; 
    doubletime=automatic, sampletime=automatic, minvalue=DEFAULT_SEEDMATRIX_MINVALUE,
)
    return _expectedseedcases(observedcases, n_seeds, doubletime, sampletime, minvalue)
end

function _expectedseedcases(
    observedcases, n_seeds, inputdoubletime, inputsampletime, minvalue
)
    doubletime = _expectedseedcasesdoubletime(inputdoubletime, n_seeds)
    sampletime = _expectedseedcasessampletime(inputsampletime, observedcases, n_seeds)
    return __expectedseedcases(observedcases, n_seeds, doubletime, sampletime, minvalue)
end

function __expectedseedcases(observedcases, n_seeds, doubletime, sampletime, minvalue)
    exptdseedcases = zeros(n_seeds, _ngroups(observedcases))

    for g in 1:_ngroups(observedcases)
        tot = sum(@view observedcases[1:sampletime, g]) / sampletime
        for t in 1:n_seeds
            exptdseedcases[t, g] = max(
                0,
                log(tot) - (n_seeds + 1 - t) * log(2) / doubletime
            )
        end
    end

    _expectedseedcasesminimum!(exptdseedcases, observedcases, minvalue)
    return exptdseedcases
end

_expectedseedcasesdoubletime(inputdoubletime, ::Any) = inputdoubletime
_expectedseedcasesdoubletime(::Automatic, n_seeds) = n_seeds
_expectedseedcasessampletime(inputsampletime, ::Any, ::Any) = inputsampletime

function _expectedseedcasessampletime(::Automatic, observedcases, n_seeds)
    return min(n_seeds, _ntimes(observedcases))
end

function _expectedseedcasesminimum!(exptdseedcases, observedcases, minvalue)
    for g in 1:_ngroups(observedcases)
        sum(exptdseedcases[:, g]) >= minvalue && continue 
        exptdseedcases[1, g] += minvalue - sum(exptdseedcases[:, g])
    end
    return nothing
end

_expectedseedcasesminimum!(::Any, ::Any, ::Nothing) = nothing

_expectedseedcasesifneeded(exptdseedcases::Matrix, ::Any, ::Any; kwargs...) = exptdseedcases

function _expectedseedcasesifneeded(::Nothing, observedcases, n_seeds; kwargs...)
    return expectedseedcases(observedcases, n_seeds; kwargs...)
end

function _expectedinfections(g, args...; kwargs...) 
    return _expectedinfections(generationtime, args...; func=g, kwargs...)
    # `_expectedinfections` will accept a function or a vector from the keyword argument `func`
end

function _expectedinfections(g::_Useablegenerationfunctions, logR_0, hx; kwargs...)
    return _expectedinfections(g, 1, logR_0, hx; kwargs...)
end

function _expectedinfections(
    g::_Useablegenerationfunctions, inputpropsus::T, logR_0, hx; 
    kwargs...
) where T
    # rather than throw an error, enforce 0 <= propsus <= 1 
    propsus = min(one(T), max(zero(T), inputpropsus))
    t = length(hx) + 1
    return sum(exp(logR_0) .* propsus .* [hx[x] * g(t - x; kwargs...) for x in eachindex(hx)])
end

_approxcasescalc(x, sigma) = x * (1 + sigma)
_approxcases(x, sigma) = __approxcases(_approxcasescalc(x, sigma))
_approxcases(x, sigma, ceiling) = __approxcases(_approxcasescalc(x, sigma), ceiling)

function __approxcases(x_oneplussigma::T) where T
    return max(zero(T), x_oneplussigma)
end

function __approxcases(x_oneplussigma::T, ceiling) where T
    return min(max(zero(T), x_oneplussigma), ceiling)
end

function _infections(
    g, M_x, logR_0::Matrix{T}, exptdseedcases, Ns, n_seeds; 
    kwargs...
) where T
    return _infections(g, T, M_x, logR_0, exptdseedcases, Ns, n_seeds; kwargs...)
end

function _infections(g, T::DataType, M_x, logR_0, exptdseedcases, Ns, n_seeds; kwargs...) 
    infn = _infectionsmatrix(T, logR_0, n_seeds)
    _infections!(g, infn, M_x, logR_0, exptdseedcases, Ns, n_seeds; kwargs...) 
    return infn
end

# `_infectionsmatrix` will always return a `Matrix{<:Complex}` and will give a warning if 
# something different is requested

function _infectionsmatrix(logR_0, n_seeds) 
    return _infectionsmatrix(ComplexF64, logR_0, n_seeds)
end

function _infectionsmatrix(T::DataType, logR_0, n_seeds)
    if T <: Complex
        return zeros(T, _ntimes(logR_0) + n_seeds, _ngroups(logR_0))
    else 
    _realinfectionmatrixwarning(T)
        return _infectionsmatrix(ComplexF64, logR_0, n_seeds)
    end
end

function _infections!(g, infn, M_x, logR_0, exptdseedcases, Ns, n_seeds; kwargs...) 
    _infectionsassertions(infn, M_x, logR_0, exptdseedcases, Ns, n_seeds)
    _infections_seed!(infn, M_x, exptdseedcases, Ns, n_seeds; kwargs...) 
    _infections_transmitted!(g, infn, M_x, logR_0, Ns, n_seeds; kwargs...)
    return nothing
end

function _infections_seed!(infn, M_x, exptdseedcases, ::Nothing, n_seeds; kwargs...) 
    return __infections_seed!(infn, M_x, exptdseedcases, n_seeds; kwargs...) 
end

# version that tracks proportion susceptible
function _infections_seed!(infn, M_x, exptdseedcases, Ns, n_seeds; kwargs...) 
    minimum(Ns) > 0 || throw(_zeropopulationerror(Ns))
    __infections_seed!(infn, M_x, exptdseedcases, n_seeds; kwargs...) 
    for j in 1:_ngroups(M_x), t in 1:n_seeds
        # previous proportion susceptible
        prevsus = _prevpropsus(infn, t, j; kwargs...) * Ns[j]  
        if real(infn[t, j]) > prevsus  # use this as a maximum 
            infn[t, j] = prevsus + 0im 
        else  # add new proportion susceptible
            infn[t, j] += ((prevsus - real(infn[t, j])) / Ns[j])im
        end
    end
    return nothing
end

function __infections_seed!(infn, M_x, exptdseedcases, n_seeds; kwargs...)
    for j in 1:_ngroups(M_x), t in 1:n_seeds
        infn[t, j] = _approxcases(exptdseedcases[t, j], M_x[t, j])
    end 
    return nothing
end

function _infections_transmitted!(g, infn, M_x, logR_0, ::Nothing, n_seeds; kwargs...) 
    for j in 1:_ngroups(M_x), t in 1:_ntimes(logR_0)
        infn[t+n_seeds, j] = _approxcases(
            _expectedinfections(
                g, logR_0[t, j], real.(@view infn[1:t+n_seeds-1, j]); 
                kwargs...
            ), 
            M_x[t+n_seeds, j]
        )
    end
    return nothing
end

function _infections_transmitted!(g, infn, M_x, logR_0, Ns, n_seeds; kwargs...) 
    for j in 1:_ngroups(M_x), t in 1:_ntimes(logR_0)
        prevsus = _prevpropsus(infn, t + n_seeds, j; kwargs...) * Ns[j] 
        expectednewcases = _approxcases(
            _expectedinfections(
                g, 
                prevsus / Ns[j], 
                logR_0[t, j], 
                real.(@view infn[1:t+n_seeds-1, j]); 
                kwargs...
            ), 
            M_x[t+n_seeds, j],
            prevsus
        )
        newcases = min(expectednewcases, prevsus)
        infn[t+n_seeds, j] = newcases + ((prevsus - newcases) / Ns[j])im
    end
    return nothing
end

function _prevpropsus(infn, t, j; kwargs...)
    if t == 1 
        return one(_prevpropsus(infn, 2, j))
    end
    propsus = imag(infn[t-1, j])
    return propsus
end

function renewaldid(
    data::AbstractRenewalDiDData, g, priors::RenewalDiDPriors{Q, R, S, T, U, V, W, X, Y}; 
    kwargs...
) where {
    Q <: Distribution, 
    R <: Distribution, 
    S <: Number, 
    T <: Number, 
    U <: Distribution, 
    V <: Distribution, 
    W <: Distribution, 
    X <: Distribution, 
    Y <: Any
}
    n_seeds = size(data.exptdseedcases, 1)
    return _renewaldid(
        data.observedcases,
        data.interventions,
        data.exptdseedcases,
        data.Ns,
        g,    
        priors.alphaprior,
        priors.M_xprior,
        priors.mu_delayprior,
        priors.psiprior,
        priors.sigma_delayprior,
        priors.sigma_gammaprior,
        priors.sigma_thetaprior,
        priors.tauprior,
        n_seeds,
        priors.omegaprior;
        kwargs...
    )
end

@model function _renewaldid(
    observedcases,
    interventions,
    expectedseedcases,
    Ns,
    g,    
    alphaprior,
    M_xprior,
    mu_delayprior,
    psiprior,
    sigma_delayprior,
    sigma_gammaprior,
    sigma_thetaprior,
    tauprior,
    n_seeds,
    omega;
    kwargs...
)
    ngroups = _ngroups(interventions)
    ntimes = _ntimes(interventions)

    tau ~ tauprior
    alpha ~ alphaprior
    sigma_gamma ~ sigma_gammaprior
    gammas_raw ~ filldist(Normal(0, 1), ngroups - 1)
    thetas_raw ~ filldist(Normal(0, 1), ntimes - 1)
    sigma_theta ~ sigma_thetaprior
    psi ~ psiprior
    #mu_delay ~ mu_delayprior
    #sigma_delay ~ sigma_delayprior
    M_x ~ filldist(truncated(Normal(0, 1); lower=-1), ntimes + n_seeds, ngroups)
    fittingsigma ~ Exponential(1)
    predictobservedinfectionssigmamatrix ~ filldist(
        truncated(Normal(0, 1); lower=-1), ntimes + 1, ngroups
    )

    # attempting to fit `mu_delay` and `sigma_delay` gives NaN gradients so currently use
    # constants
    mu_delay = log(2)
    sigma_delay = log(5)

    _z = max(
        maximum(
            [
                [
                    tau, 
                    alpha, 
                    sigma_gamma, 
                    sigma_theta, 
                    psi, 
                    mu_delay, 
                    sigma_delay, 
                    fittingsigma
                ]; 
                gammas_raw; 
                thetas_raw
            ]
        ),
        maximum(M_x),
        maximum(predictobservedinfectionssigmamatrix)
    ) 

    gammavec = _gammavec(gammas_raw, sigma_gamma)
    thetavec = _thetavec(thetas_raw, sigma_theta)

    predictedlogR_0 = _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)

    T = Complex{typeof(predictedlogR_0[1, 1])}
    predictedinfections = _infectionsmatrix(T, predictedlogR_0, n_seeds)
    _infections!(
        g, predictedinfections, M_x, predictedlogR_0, expectedseedcases, Ns, n_seeds; 
        kwargs...
    )

    # equal delay assumed for all cases 
    delaydistn = LogNormal(mu_delay, sigma_delay)
    delayedinfections = zeros(T, n_seeds + ntimes, ngroups)

    for t in 1:(n_seeds + ntimes) 
        for j in 1:ngroups
            newvalue = zero(T) 
            for x in 1:t 
                newvalue += _delayedinfections(predictedinfections, delaydistn, x, t, j)
            end
            delayedinfections[t, j] += newvalue 
        end
    end

    # Normal approximation of Binomial to avoid forcing integer values 
    np = real.(delayedinfections[n_seeds:n_seeds+ntimes, :]) .* psi

    if isnan(maximum(np)) 
        @addlogprob! -Inf
        return  # exit the model evaluation early
    end

    predictobservedinfections = max.(0, np .+ predictobservedinfectionssigmamatrix .* np .* (1 - psi))

    observedcases ~ arraydist(Normal.(predictobservedinfections, fittingsigma))
end

function _delayedinfections(
    predictedinfections::AbstractArray{T}, delaydistn, x, t, j
) where T
    if isnan(cdf(delaydistn, t + 1 - x)) || isnan(cdf(delaydistn, t - x))
        return zero(T)
    end
    return predictedinfections[x, j] * (cdf(delaydistn, t + 1 - x) - cdf(delaydistn, t - x))
end


# Assertions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _diddataassertions(observedcases, interventions, exptdseedcases) 
    _ngroups(observedcases) == _ngroups(interventions) || throw(
        _ngroupsdimensionmismatch("observedcases", "interventions")
    )
    _ngroups(observedcases) == _ngroups(exptdseedcases) || throw(
        _ngroupsdimensionmismatch("observedcases", "exptdseedcases")
    )
    _ntimes(observedcases) == _ntimes(interventions) + 1 || throw(_ntimesdimensionmismatch())
    return nothing
end

function _diddataassertions(observedcases, interventions, Ns, exptdseedcases) 
    _diddataassertions(observedcases, interventions, exptdseedcases) 
    _ngroups(observedcases) == length(Ns) || throw(
        _ngroupsdimensionmismatch("observedcases", "Ns")
    )
    return nothing
end

function _infectionsassertions(infn, M_x, logR_0, exptdseedcases, ::Nothing, n_seeds)
    _expctdtotaltimes = _ntimes(logR_0) + _ntimes(exptdseedcases)
    _ngroups(infn) == _ngroups(logR_0) || throw(_widthmismatch("infn", "logR_0"))
    _ntimes(infn) == _expctdtotaltimes || throw(_infectionsntimeserror("infn"))
    _ngroups(M_x) == _ngroups(logR_0) || throw(_widthmismatch("M_x", "logR_0"))
    _ngroups(M_x) == _ngroups(exptdseedcases) || throw(_widthmismatch("M_x", "exptdseedcases"))
    _ntimes(M_x) == _expctdtotaltimes || throw(_infectionsntimeserror("M_x"))
    _ntimes(exptdseedcases) == n_seeds || throw(_infectionsnseedserror())
    return nothing
end

function _infectionsassertions(infn, M_x, logR_0, exptdseedcases, Ns, n_seeds)
    _infectionsassertions(infn, M_x, logR_0, exptdseedcases, nothing, n_seeds)
    _ngroups(M_x) == length(Ns) || throw(_MxNserror(Ns, M_x))
    return nothing
end

function _predictedlogR_0assertions(gammas, thetas, interventions)
    length(gammas) == _ngroups(interventions) || throw(_ngroupsmmerror(gammas, interventions))
    length(thetas) == _ntimes(interventions) || throw(_ntimesmmerror(thetas, interventions))
    return nothing
end


# Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _realinfectionmatrixwarning(T)
    m = "`infectionsmatrix` requested with type `Matrix{$T}`. Must be complex so type \
        `Matrix{ComplexF64}` returned"
    @warn m
    return nothing 
end

# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _infectionsnseedserror()
    return DimensionMismatch("height of `exptdseedcases` must equal `n_seeds`")
end

function _infectionsntimeserror(M)
    m = "height of `$M` must equal sum of heights of `logR_0` and `exptdseedcases`"
    return DimensionMismatch(m)
end

_MxNserror(Ns, M_x) = _vm_dimensionmismatch(length(Ns), "Ns", "width", "M_x", _ngroups(M_x))

function _ngroupsdimensionmismatch(a, b)
    return DimensionMismatch("number of groups must be equal for $a and $b")
end

function _ngroupsmmerror(gammavec, interventions)
    return _R_0_dimensionmismatch(
        length(gammavec), "gammavec", "width", _ngroups(interventions)
    )
end

function _ntimesdimensionmismatch()
    m = "number of times of observed cases must be 1 greater than number of times for \
        interventions"
    return DimensionMismatch(m)
end

function _ntimesmmerror(thetavec, interventions)
    return _R_0_dimensionmismatch(
        length(thetavec), "thetavec", "height", _ntimes(interventions)
    )
end

_propsusargumenterror(propsus) = ArgumentError("$propsus: must satisfy 0 ≤ propsus ≤ 1")

function _R_0_dimensionmismatch(len, vec, dim, size)
    return _vm_dimensionmismatch(len, vec, dim, "interventions", size)
end

function _vm_dimensionmismatch(len, vec, dim, mat, size)
    m = "length of `$vec`, $len, should equal $dim of `$mat`, $size"
    return DimensionMismatch(m)
end

_widthmismatch(a, b) = DimensionMismatch("`$a` and `$b` must have the same width")

function _zeropopulationerror(Ns)
    m = "$Ns: all population sizes must be greater than 0. To use the function without a \
        population size use `Ns=nothing`"
    return ArgumentError(m)
end
