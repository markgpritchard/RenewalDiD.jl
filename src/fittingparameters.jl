# Use Turing.jl to fit parameters to a dataset

## Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

## Functions called by `_renewaldid`

_ntimes(M::AbstractMatrix) = size(M, 1)
_ngroups(M::AbstractMatrix) = size(M, 2)

# the mean of all gamma values is zero, ensured by setting the final value in the vector as 
# `-sum` of other values
_gammavec(gammas_raw, sigma_gamma) = [gammas_raw; -sum(gammas_raw)] .* sigma_gamma

# The value of theta at time 0 is fixed to 0. `thetas_raw` is a vector representing how much 
# each subsequent theta differs from the previous one as a multiple of the standard 
# deviation. `sigma_theta` is the standard deviation. The output of this function is the 
# time-varying `theta` at each time as a cumulative random walk. 
function _thetavec(thetas_raw::AbstractVector{S}, sigma_theta::T) where {S, T}
    theta_0 = zero(S) * zero(T)
    return cumsum([theta_0; thetas_raw .* sigma_theta])
end

function _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)
    _predictedlogR_0assertions(gammavec, thetavec, interventions)

    R_0 = alpha .* ones(size(interventions)...)  # intercept
    
    for (g, gamma) in enumerate(gammavec)  # group-dependent values 
        R_0[:, g] .+=  gamma
    end 
    
    for (t, theta) in enumerate(thetavec)  # time-dependent values 
        R_0[t, :] .+=  theta
    end
    
    R_0 .+= tau .* interventions  # effect of the intervention
    return R_0
end

function _expectedseedcases(observedcases, n_seeds; doubletime=n_seeds, sampletime=automatic)
    return __expectedseedcases(observedcases, n_seeds, doubletime, sampletime)
end

function __expectedseedcases(observedcases, n_seeds, doubletime, ::Automatic)
    sampletime = min(n_seeds, _ntimes(observedcases))
    return __expectedseedcases(observedcases, n_seeds, doubletime, sampletime)
end

function __expectedseedcases(observedcases, n_seeds, doubletime, sampletime)
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

    return exptdseedcases
end

function _expectedinfections(g, logR_0, hx; kwargs...)
    t = length(hx) + 1
    return sum(exp(logR_0) .* [hx[x] * g(t - x; kwargs...) for x in eachindex(hx)])
end

_approxcases(x, sigma) = max(0, x * (1 + sigma))

function _infections(
    g, M_x, logR_0::Matrix{T}, exptdseedcases, Ns, n_seeds; 
    kwargs...
) where T
    _infectionsassertions(M_x, logR_0, exptdseedcases, Ns, n_seeds)
    infn = zeros(T, _ntimes(logR_0) + n_seeds, _ngroups(logR_0))

    for j in 1:_ngroups(M_x), t in 1:n_seeds
        infn[t, j] = _approxcases(exptdseedcases[t, j], M_x[t, j])
    end
    
    for j in 1:_ngroups(M_x), t in 1:_ntimes(logR_0)
        infn[t+n_seeds, j] = _approxcases(
            _expectedinfections(g, logR_0[t, j], @view infn[1:t+n_seeds-1, j]; kwargs...), 
            M_x[t+n_seeds, j]
        )
    end
    
    return infn
end

function packdata(; observedcases, interventions, Ns)
    return Dict(
        :observedcases => observedcases,
        :interventions => interventions,
        :Ns => Ns,
    )
end

function packpriors( ;
    alphaprior=Normal(0, 1),
    mu_gammaprior=Normal(0, 1),
    sigma_gammaprior=Exponential(1),
    sigma_thetaprior=Exponential(1),
    tauprior=Normal(0, 1),
)
    return Dict(
        :alphaprior => alphaprior,
        :mu_gammaprior => mu_gammaprior,
        :sigma_gammaprior => sigma_gammaprior,
        :sigma_thetaprior => sigma_thetaprior,
        :tauprior => tauprior,
    )
end

function renewaldid(
    data, g, priors; 
    n_seeds=7, doubletime=n_seeds, sampletime=automatic, kwargs...
)
    @unpack observedcases, interventions, Ns = data 
    @unpack alphaprior, mu_gammaprior, sigma_gammaprior, sigma_thetaprior, tauprior = priors
    return _renewaldid(
        observedcases,
        interventions,
        Ns,
        g,    
        alphaprior,
        mu_gammaprior,
        sigma_gammaprior,
        sigma_thetaprior,
        tauprior,
        n_seeds,
        doubletime, 
        sampletime;
        kwargs...
    )
end

@model function _renewaldid(
    observedcases,
    interventions,
    Ns,
    g,    
    alphaprior,
    mu_gammaprior,
    sigma_gammaprior,
    sigma_thetaprior,
    tauprior,
    n_seeds,
    doubletime, 
    sampletime;
    kwargs...
)
    ngroups = _ngroups(interventions)
    ntimes = _ntimes(interventions)

    tau ~ tauprior
    alpha ~ alphaprior
    mu_gamma ~ mu_gammaprior
    sigma_gamma ~ sigma_gammaprior
    gammas_raw ~ filldist(Normal(0, 1), ngroups - 1)
    thetas_raw ~ filldist(Normal(0, 1), ntimes - 1)
    sigma_theta ~ sigma_thetaprior
    M_x ~ filldist(Normal(0, 1), ntimes + n_seeds, ngroups)

    gammavec = _gammavec(gammas_raw, sigma_gamma)
    thetavec = _thetavec(thetas_raw, sigma_theta)

    predictedlogR_0 = _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)

    expectedseedcases = _expectedseedcases(observedcases, n_seeds; doubletime, sampletime)
    predictedinfections = _infections(
        g, M_x, predictedlogR_0, expectedseedcases, Ns, n_seeds; 
        kwargs...
    )

    # to add delay later 
    observedcases ~ arraydist(Normal.(predictedinfections[n_seeds:n_seeds+ntimes, :], 1))
end



# Assertions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _infectionsassertions(M_x, logR_0, exptdseedcases, Ns, n_seeds)
    _ngroups(M_x) == _ngroups(logR_0) || throw(_widthmismatch("M_x", "logR_0"))
    _ngroups(M_x) == _ngroups(exptdseedcases) || throw(_widthmismatch("M_x", "exptdseedcases"))
    _ngroups(M_x) == length(Ns) || throw(_MxNserror(Ns, M_x))
    _ntimes(M_x) == _ntimes(logR_0) + _ntimes(exptdseedcases) || throw(_infectionsntimeserror())
    _ntimes(exptdseedcases) == n_seeds || throw(_infectionsnseedserror())
    return nothing
end

function _predictedlogR_0assertions(gammas, thetas, interventions)
    length(gammas) == _ngroups(interventions) || throw(_ngroupsmmerror(gammas, interventions))
    length(thetas) == _ntimes(interventions) || throw(_ntimesmmerror(thetas, interventions))
    return nothing
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _infectionsnseedserror()
    return DimensionMismatch("height of `exptdseedcases` must equal `n_seeds`")
end

function _infectionsntimeserror()
    return DimensionMismatch(
        "height of `Mx` must equal sum of heights of `logR_0` and `exptdseedcases`"
    )
end

_MxNserror(Ns, M_x) = _vm_dimensionmismatch(length(Ns), "Ns", "width", "M_x", _ngroups(M_x))

function _ngroupsmmerror(gammavec, interventions)
    return _R_0_dimensionmismatch(
        length(gammavec), "gammavec", "width", _ngroups(interventions)
    )
end

function _ntimesmmerror(thetavec, interventions)
    return _R_0_dimensionmismatch(
        length(thetavec), "thetavec", "height", _ntimes(interventions)
    )
end

function _R_0_dimensionmismatch(len, vec, dim, size)
    return _vm_dimensionmismatch(len, vec, dim, "interventions", size)
end

function _vm_dimensionmismatch(len, vec, dim, mat, size)
    return DimensionMismatch("length of `$vec`, $len, should equal $dim of `$mat`, $size")
end

_widthmismatch(a, b) = DimensionMismatch("`$a` and `$b` must have the same width")
