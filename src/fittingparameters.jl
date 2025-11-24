# Use Turing.jl to fit parameters to a dataset

# Structs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## struct of prior distributions 

"""
    RenewalDiDPriors{Q, S, T, U, V, W, X}

Container for prior distributions passed to `renewaldid`.

# Fields
- `alphaprior::Q=Normal(0, 1)`: intercept for log basic reproduction ratio 
- `sigma_gammaprior::S=Exponential(1)`: variance in group-specific contributions to log 
    basic reproduction ratio
- `sigma_thetaprior::T=Exponential(1)`: variance in time-varying contributions to log 
    basic reproduction ratio
- `tauprior::U=Normal(0, 1)`: average treatment effect in the treated
- `psiprior::V=Beta(1, 1)`: proportion of infections that are diagnosed
- `omegaprior::W=0`: rate of return to susceptibility after infection (not currently used)
- `delaydistn::X=LogNormal(log(2), log(2))`: distribution of delay times between infection 
    and diagnosis

All arguments are optional with default values as above. The constructor takes each prior as 
    a keyword argument. 

# Examples
```jldoctest
julia> using Turing

julia> RenewalDiDPriors(; psiprior=Beta(6, 4))
RenewalDiDPriors{Normal{Float64}, Exponential{Float64}, Exponential{Float64}, \
    Normal{Float64}, Beta{Float64}, Int64, LogNormal{Float64}}
 alphaprior:       Normal{Float64}(μ=0.0, σ=1.0)
 sigma_gammaprior: Exponential{Float64}(θ=1.0)
 sigma_thetaprior: Exponential{Float64}(θ=1.0)
 tauprior:         Normal{Float64}(μ=0.0, σ=1.0)
 psiprior:         Beta{Float64}(α=6.0, β=4.0)
 omegaprior:       0
 delaydistn:       LogNormal{Float64}(μ=0.6931471805599453, σ=0.6931471805599453)
```
""" 
@kwdef struct RenewalDiDPriors{Q, S, T, U, V, W, X}
    # attempting to fit a delay distribution has so far given NaN gradients. Function 
    # `renewaldid` currently expects a distribution to be supplied 
    alphaprior::Q=Normal(0, 1)
    sigma_gammaprior::S=Exponential(1)
    sigma_thetaprior::T=Exponential(1)
    tauprior::U=Normal(0, 1)
    psiprior::V=Beta(1, 1)
    omegaprior::W=0
    delaydistn::X=LogNormal(log(2), log(2))
end

function Base.show(
    io::IO, ::MIME"text/plain", p::RenewalDiDPriors{Q, S, T, U, V, W, X}
) where {Q, S, T, U, V, W, X}
    print(
        io,
        "RenewalDiDPriors{$Q, $S, $T, $U, $V, $W, $X}",
        "\n alphaprior:       ", p.alphaprior,
        "\n sigma_gammaprior: ", p.sigma_gammaprior,
        "\n sigma_thetaprior: ", p.sigma_thetaprior,
        "\n tauprior:         ", p.tauprior,
        "\n psiprior:         ", p.psiprior,
        "\n omegaprior:       ", p.omegaprior,
        "\n delaydistn:       ", p.delaydistn
    )
end


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

## Functions called by `_renewaldid`

_ntimes(A::AbstractArray) = size(A, 1)
_ngroups(A::AbstractArray) = size(A, 2)
_ninterventions(A::AbstractArray) = size(A, 3)
_nthetas(ntimes, ::Automatic) = ntimes - 1
_nthetas(ntimes, thetainterval::Integer) = round(Int, (ntimes - 1) / thetainterval, RoundUp)

# the mean of all gamma values is zero, ensured by setting the final value in the vector as 
# `-sum` of other values
_gammavec(gammas_raw, sigma_gamma) = [gammas_raw; -sum(gammas_raw)] .* sigma_gamma
# a second version of this function is given in `processparameters.jl`

# The value of theta at time 0 is fixed to 0. `thetas_raw` is a vector representing how much 
# each subsequent theta differs from the previous one as a multiple of the standard 
# deviation. `sigma_theta` is the standard deviation. The output of this function is the 
# time-varying `theta` at each time as a cumulative random walk. 
function _thetavec(thetas_raw::AbstractVector, sigma_theta, ntimes; thetainterval=automatic) 
    rawthetavec = _assemblethetavec(thetas_raw, sigma_theta, ntimes, thetainterval) 
    return cumsum(rawthetavec)
end  # a second version of this function is given in `processparameters.jl`

function _assemblethetavec(
    thetas_raw, 
    sigma_theta, 
    ntimes, 
    ::Automatic, 
    theta_0=_thetaveczero(thetas_raw, sigma_theta)
)
    length(thetas_raw) == ntimes - 1 || throw(_thetaveclengthmismatch(thetas_raw, ntimes))
    return [theta_0; thetas_raw .* sigma_theta]
end

function _assemblethetavec(
    thetas_raw, 
    sigma_theta, 
    ntimes, 
    thetainterval::Integer, 
    theta_0=_thetaveczero(thetas_raw, sigma_theta, thetainterval)
) 
    if thetainterval == 1  
        return _assemblethetavec(thetas_raw, sigma_theta, ntimes, automatic, theta_0)
    elseif thetainterval <= 0 
        throw(_assemblethetavecnonintegererror(thetainterval))
    else
        return __assemblethetavec(thetas_raw, sigma_theta, ntimes, thetainterval, theta_0)
    end
end

function __assemblethetavec(thetas_raw, sigma_theta, ntimes, thetainterval, ::T) where T
    _assemblethetavecassertion(thetas_raw, ntimes, thetainterval)
    thetavec = zeros(T, ntimes)
    k = 1 
    ℓ = zero(T)

    for (j, v) in enumerate(thetas_raw)
        k >= ntimes && continue  # avoid trying to allocate beyond end of vector
        thetavec[k] = ℓ * sigma_theta
        k += 1
        for n in 1:(thetainterval - 1)
            k >= ntimes && continue  # avoid trying to allocate beyond end of vector
            thetavec[k] = (ℓ * (thetainterval - n) + v * n) * sigma_theta / thetainterval 
            k += 1
        end
        ℓ = v 
    end

    thetavec[ntimes] = last(thetas_raw) * sigma_theta
    return thetavec
end

_thetaveczero(::AbstractVector{S}, ::T) where {S, T} = zero(S) * zero(T)
_thetaveczero(::AbstractVector{S}, ::T, ::U) where {S, T, U} = zero(S) * zero(T) / one(U)

function _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)
    _predictedlogR_0assertions(gammavec, thetavec, interventions)
    logR_0 = [
        alpha + gammavec[g] + thetavec[t] + _totaltau(tau, interventions, t, g) 
        for t in axes(interventions, 1), g in axes(interventions, 2)
    ]
    return logR_0
end

_totaltau(tau, interventions::AbstractMatrix, t, g) = tau[1] * interventions[t, g]

function _totaltau(tau, interventions::AbstractArray{<:Any, 3}, t, g) 
    return sum([tau[k] * interventions[t, g, k] for k in axes(interventions, 3)])
end

"""
    expectedseedcases(observedcases, n_seeds=7; <keyword arguments>)

Return a matrix of infections before time `0` to seed subsequent transmission.

# Arguments 
- `observedcases` is the matrix of observations that is used in the analysis.
- `n_seeds=7` is the number of rows in the matrix of seed cases.

# Keyword arguments 
- `doubletime`: assumed doubling time for infections before time `0`; the default is to 
    equal `n_seeds`
- `sampletime`: how many rows of `observedcases` to consider when deciding on the number of 
    infections to include in the seed matrix; default is `doubletime` if provided, otherwise 
    `n_seeds`
- `minvalue=0.5`: minimum number of infection to include in seed matrix regardless of 
    numbers in `observedcases`

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(10);

julia> observedcases = rand(rng, Poisson(10), 11, 3);

julia> expectedseedcases(observedcases, 3)
3×3 Matrix{Float64}:
 1.7346   1.42712  1.79176
 1.96565  1.65817  2.02281
 2.1967   1.88921  2.25386

julia> expectedseedcases(observedcases, 3; doubletime=5)
3×3 Matrix{Float64}:
 2.01186  1.70438  2.06902
 2.15049  1.843    2.20765
 2.28912  1.98163  2.34628

julia> expectedseedcases(observedcases, 3; minvalue=10)
3×3 Matrix{Float64}:
 1.7346   1.42712  1.79176
 1.96565  1.65817  2.02281
 6.29975  6.91472  6.18543
```
"""
function expectedseedcases(
    observedcases, n_seeds=7; 
    doubletime=automatic, sampletime=automatic, minvalue=0.5,
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
    _expectedseedcases!(exptdseedcases, observedcases, n_seeds, doubletime, sampletime)
    _expectedseedcasesminimum!(exptdseedcases, observedcases, minvalue, n_seeds)
    return exptdseedcases
end

function _expectedseedcases!(exptdseedcases, observedcases, n_seeds, doubletime, sampletime)
    for g in 1:_ngroups(observedcases)
        tot = sum(@view observedcases[1:sampletime, g]) / sampletime
        for t in 1:n_seeds
            exptdseedcases[t, g] = exp(log(tot) - (n_seeds + 1 - t) * log(2) / doubletime)
        end
    end
    return nothing
end

_expectedseedcasesdoubletime(inputdoubletime, ::Any) = inputdoubletime
_expectedseedcasesdoubletime(::Automatic, n_seeds) = n_seeds
_expectedseedcasessampletime(inputsampletime, ::Any, ::Any) = inputsampletime

function _expectedseedcasessampletime(::Automatic, observedcases, n_seeds)
    return min(n_seeds, _ntimes(observedcases))
end

function _expectedseedcasesminimum!(exptdseedcases, observedcases, minvalue, n_seeds)
    for g in 1:_ngroups(observedcases)
        c = sum(@view exptdseedcases[:, g])
        c >= minvalue && continue 
        exptdseedcases[n_seeds, g] += minvalue - c
    end
    return nothing
end

_expectedseedcasesminimum!(::Any, ::Any, ::Nothing, ::Any) = nothing
_expectedseedcasesifneeded(exptdseedcases::Matrix, ::Any; kwargs...) = exptdseedcases

function _expectedseedcasesifneeded(::Nothing, observedcases; n_seeds=automatic, kwargs...)
    return __expectedseedcasesifneeded(observedcases, n_seeds; kwargs...)
end

function __expectedseedcasesifneeded(observedcases, ::Automatic; kwargs...)
    return expectedseedcases(observedcases; kwargs...)
end

function __expectedseedcasesifneeded(observedcases, n_seeds; kwargs...)
    return expectedseedcases(observedcases, n_seeds; kwargs...)
end

function _expectedinfections(g, args...; kwargs...) 
    # `_expectedinfections` will accept a function or a vector in the keyword argument `func`
    return _expectedinfections(generationtime, args...; func=g, kwargs...)
end

function _expectedinfections(g::_Useablegenerationfunctions, logR_0, hx; kwargs...)
    return __expectedinfections(g, 1, logR_0, hx; kwargs...)  # assume propsus = 1
end

function _expectedinfections(
    g::_Useablegenerationfunctions, inputpropsus, logR_0, hx; 
    kwargs...
) 
    propsus = min(1, max(0, inputpropsus))  # enforce 0 <= propsus <= 1 
    return __expectedinfections(g, propsus, logR_0, hx; kwargs...) 
end

function __expectedinfections(g, propsus, logR_0, hx; kwargs...) 
    t = length(hx) + 1
    return sum(exp(logR_0) .* propsus .* [hx[x] * g(t - x; kwargs...) for x in eachindex(hx)])
end

_approxcasescalc(x, sigma) = __approxcasescalc(max(zero(x), x), sigma)  # should never x < 0
__approxcasescalc(x, sigma) = x + sigma * NaNMath.sqrt(x)

# put a ceiling in all instances to prevent NaN outcomes when parameters generate
# astronomical numbers of infections
function _approxcases(x, sigma, ceiling=1e12)
    calculatedcases = _approxcasescalc(x, sigma)
    # if isnan(calculated), return ceiling
    return NaNMath.min(max(zero(calculatedcases), calculatedcases), ceiling)
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

_infectionsmatrix(logR_0, n_seeds) = _infectionsmatrix(ComplexF64, logR_0, n_seeds)

function _infectionsmatrix(T::Type{<:Complex}, logR_0, n_seeds)
    return zeros(T, _ntimes(logR_0) + n_seeds, _ngroups(logR_0))
end

function _infectionsmatrix(T::DataType, logR_0, n_seeds)
    _realinfectionmatrixwarning(T)
    return _infectionsmatrix(ComplexF64, logR_0, n_seeds)
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
    return infn
end

function __infections_seed!(infn, M_x, exptdseedcases, n_seeds; kwargs...)
    for j in 1:_ngroups(M_x), t in 1:n_seeds
        infn[t, j] = _approxcases(exptdseedcases[t, j], M_x[t, j])
    end 
    return infn
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

"""
    renewaldid(data::AbstractRenewalDiDData, g, priors::RenewalDiDPriors; [thetainterval], <keyword arguments>)

Return a `Turing.@model` for parameter fitting.

`g` is a generation interval function and the keyword arguments are passed to `g`.

`thetainterval` is the interval between independently sampled time-varying parameters. If 
    `thetainterval==1` (the default) a new value is entered each day.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(1);

julia> sim = RenewalDiD.testsimulation(rng);

julia> renewaldid(sim, g_seir, RenewalDiDPriors(); mu=0.2, kappa=0.5)
DynamicPPL.Model{typeof(RenewalDiD._renewaldid), (:observedcases, :interventions, \
    :expectedseedcases, :Ns, :g, :alphaprior, :psiprior, :sigma_gammaprior, \
    :sigma_thetaprior, :tauprior, :delaydistn, :n_seeds, :omega), (:mu, :kappa), (), \
    Tuple{Matrix{Int64}, InterventionMatrix{Int64}, Matrix{Float64}, Vector{Int64}, \
    typeof(g_seir), Normal{Float64}, Beta{Float64}, Exponential{Float64}, \
    Exponential{Float64}, Normal{Float64}, LogNormal{Float64}, Int64, Int64}, \
    Tuple{Float64, Float64}, DynamicPPL.DefaultContext}(RenewalDiD._renewaldid, \
    (observedcases = [0 0 0; 0 0 0; … ; 0 4 1; 1 1 4], interventions = [0 0 0; 0 0 0; … ; \
    0 1 1; 0 1 1] {duration 10, starttimes [nothing, 4, 6]}, expectedseedcases = [0.0 0.0 \
    0.0; 0.0 0.0 0.0; … ; 0.0 0.0 0.0; 0.5 0.5 0.5], Ns = [100, 200, 50], g = g_seir, \
    alphaprior = Distributions.Normal{Float64}(μ=0.0, σ=1.0), psiprior = \
    Distributions.Beta{Float64}(α=1.0, β=1.0), sigma_gammaprior = \
    Distributions.Exponential{Float64}(θ=1.0), sigma_thetaprior = \
    Distributions.Exponential{Float64}(θ=1.0), tauprior = \
    Distributions.Normal{Float64}(μ=0.0, σ=1.0), delaydistn = \
    Distributions.LogNormal{Float64}(μ=0.6931471805599453, \
    σ=0.6931471805599453), n_seeds = 7, omega = 0), (mu = 0.2, kappa = 0.5), \
    DynamicPPL.DefaultContext())
```
"""
function renewaldid(
    data::AbstractRenewalDiDData, g, priors::RenewalDiDPriors{Q, S, T, U, V, W, X}; 
    thetainterval=automatic, kwargs...
) where {
    Q <: Distribution, 
    S <: Distribution, 
    T <: Distribution, 
    U <: Distribution, 
    V <: Distribution, 
    W <: Real, 
    X <: Distribution, 
}
    return _renewaldid(
        _observedcases(data),
        _interventions(data),
        _expectedseedcases(data),
        _ns(data),
        g,    
        priors.alphaprior,
        priors.psiprior,
        priors.sigma_gammaprior,
        priors.sigma_thetaprior,
        priors.tauprior,
        priors.delaydistn,
        _nseeds(data),
        priors.omegaprior,
        thetainterval;
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
    psiprior,
    sigma_gammaprior,
    sigma_thetaprior,
    tauprior,
    delaydistn,
    n_seeds,
    omega,
    thetainterval;
    kwargs...
)
    ngroups = _ngroups(interventions)
    ntimes = _ntimes(interventions)
    ninterventions = _ninterventions(interventions)
    nthetas = _nthetas(ntimes, thetainterval)

    tau ~ filldist(tauprior, ninterventions)
    alpha ~ alphaprior
    sigma_gamma ~ sigma_gammaprior
    gammas_raw ~ filldist(Normal(0, 1), ngroups - 1)
    thetas_raw ~ filldist(Normal(0, 1), ntimes - 1)
    sigma_theta ~ sigma_thetaprior
    psi ~ psiprior
    M_x ~ filldist(Normal(0, 1), ntimes + n_seeds, ngroups)
    minsigma2 ~ Beta(1, 2)

    gammavec = _gammavec(gammas_raw, sigma_gamma)
    thetavec = _thetavec(thetas_raw, sigma_theta, ntimes; thetainterval)
    predictedlogR_0 = _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)

    T = Complex{typeof(predictedlogR_0[1, 1])}
    predictedinfections = _infectionsmatrix(T, predictedlogR_0, n_seeds)
    _infections!(
        g, predictedinfections, M_x, predictedlogR_0, expectedseedcases, Ns, n_seeds; 
        kwargs...
    )

    # delay between infection and detection
    delayedinfections = _delayedinfections(
        T, predictedinfections, delaydistn, ngroups, ntimes, n_seeds
    )

    # Normal approximation of Binomial to avoid forcing integer values 
    np = real.(delayedinfections[n_seeds:n_seeds+ntimes, :]) .* psi

    observedcases ~ arraydist(Normal.(np, NaNMath.sqrt.(np .* (1 - psi) .+ minsigma2)))
    return nothing
end

function _delayedinfections(T, predictedinfections, delaydistn, ngroups, ntimes, n_seeds)
    delayedinfections = zeros(T, n_seeds + ntimes, ngroups)
    _delayedinfections!(
        delayedinfections, predictedinfections, delaydistn, ngroups, ntimes, n_seeds
    )
    return delayedinfections
end

function _delayedinfections!(
    delayedinfections, predictedinfections, delaydistn, ngroups, ntimes, n_seeds
) 
    for t in 1:(n_seeds + ntimes) 
        for j in 1:ngroups
            newvalue = sum(
                [__delayedinfections(predictedinfections, delaydistn, x, t, j) for x in 1:t]
            )
            delayedinfections[t, j] += newvalue 
        end
    end
    return nothing
end

function __delayedinfections(predictedinfections, distribution, x, t, j) 
    return predictedinfections[x, j] * _dailyproportiondetected(t - x; distribution)
end

function _dailyproportiondetected(t::Integer; distribution)
    return cdf(distribution, t) - cdf(distribution, t - 1)
end


# Assertions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _assemblethetavecassertion(thetas_raw, ntimes, thetainterval)
    leng = length(thetas_raw) 
    requiredlength = _nthetas(ntimes, thetainterval)
    leng == requiredlength || throw(_assemblethetaveclengthrerror(leng, requiredlength))
    return nothing 
end

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

function _assemblethetaveclengthrerror(leng, requiredlength)
    return ArgumentError("$leng): length of thetas_raw must be $requiredlength")
end

function _assemblethetavecnonintegererror(thetainterval)
    return ArgumentError("$thetainterval: thetainterval must be a positive integer") 
end

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

function _thetaveclengthmismatch(thetas_raw, ntimes)
    len = length(thetas_raw) 
    expted = ntimes - 1
    return DimensionMismatch("thetas_raw length = $len, expected $expted")
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
