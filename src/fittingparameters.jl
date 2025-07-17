# Use Turing.jl to fit parameters to a dataset

## Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

## Functions called by `_renewaldid`

_ntimes(M::AbstractMatrix) = size(M, 1)
_ngroups(M::AbstractMatrix) = size(M, 2)

_gammavec(mu_gamma, sigma_gamma, gammas_raw) = mu_gamma .+ sigma_gamma .* gammas_raw

# `theta_0` is the value of theta at time 0. `thetas_raw` is a vector representing how much 
# each subsequent theta differs from the previous one as a multiple of the standard 
# deviation. `sigma_theta` is the standard deviation. The output of this function is the 
# time-varying `theta` at each time as a cumulative random walk. 
_thetavec(theta_0, thetas_raw, sigma_theta) = cumsum([theta_0; thetas_raw .* sigma_theta])

function _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)
    length(gammavec) == _ngroups(interventions) || throw(_ngroupsmmerror(gammavec, interventions))
    length(thetavec) == _ntimes(interventions) || throw(_ntimesmmerror(thetavec, interventions))
    R_0 = alpha .* ones(size(interventions)...) 
    for (g, gamma) in enumerate(gammavec)
        R_0[:, g] .+=  gamma
    end 
    for (t, theta) in enumerate(thetavec)
        R_0[t, :] .+=  theta
    end
    R_0 .+= tau .* interventions
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
    _ngroups(M_x) == _ngroups(logR_0) || throw(_widthmismatch("M_x", "logR_0"))
    _ngroups(M_x) == _ngroups(exptdseedcases) || throw(_widthmismatch("M_x", "exptdseedcases"))
    _ngroups(M_x) == length(Ns) || throw(_MxNserror(Ns, M_x))
    _ntimes(M_x) == _ntimes(logR_0) + _ntimes(exptdseedcases) || throw(_infectionsntimeserror())
    _ntimes(exptdseedcases) == n_seeds || throw(_infectionsnseedserror())
    
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
    theta_0prior=Normal(0, 1),
    sigma_thetaprior=Exponential(1),
    tauprior=Normal(0, 1),
)
    return Dict(
        :alphaprior => alphaprior,
        :mu_gammaprior => mu_gammaprior,
        :sigma_gammaprior => sigma_gammaprior,
        :theta_0prior => theta_0prior,
        :sigma_thetaprior => sigma_thetaprior,
        :tauprior => tauprior,
    )
end

function renewaldid(
    data, g, priors=packpriors(); 
    n_seeds=7, doubletime=n_seeds, sampletime=automatic, kwargs...
)
    @unpack observedcases, interventions, Ns = data 
    @unpack alphaprior, mu_gammaprior, sigma_gammaprior, theta_0prior = priors
    @unpack sigma_thetaprior, tauprior = priors
    return _renewaldid(
        observedcases,
        interventions,
        Ns,
        g,    
        alphaprior,
        mu_gammaprior,
        sigma_gammaprior,
        theta_0prior,
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
    theta_0prior,
    sigma_thetaprior,
    tauprior,
    n_seeds,
    doubletime, 
    sampletime;
    kwargs...
)
    ngroups = _ngroups(interventions)
    ntimes = _ntimes(interventions)

    alpha ~ alphaprior
    mu_gamma ~ mu_gammaprior
    sigma_gamma ~ sigma_gammaprior
    gammas_raw ~ filldist(Normal(0, 1), ngroups)
    theta_0 ~ theta_0prior
    thetas_raw ~ filldist(Normal(0, 1), ntimes - 1)
    sigma_theta ~ sigma_thetaprior
    tau ~ tauprior
    M_x ~ filldist(Normal(0, 1), ntimes + n_seeds, ngroups)

    gammavec = _gammavec(mu_gamma, sigma_gamma, gammas_raw)
    thetavec = _thetavec(theta_0, thetas_raw, sigma_theta)

    predictedlogR_0 = _predictedlogR_0(alpha, gammavec, thetavec, tau, interventions)

    expectedseedcases = _expectedseedcases(observedcases, n_seeds; doubletime, sampletime)
    predictedinfections = _infections(
        g, M_x, predictedlogR_0, expectedseedcases, Ns, n_seeds; 
        kwargs...
    )

    # to add delay later 
    observedcases ~ arraydist(Normal.(predictedinfections[n_seeds:n_seeds+ntimes, :], 1))
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
