# run simulation with fitted parameters 

## Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

function samplerenewaldidinfections(g, df, i; interventions, Ns, seedmatrix, ngroups, ntimes)
    return samplerenewaldidinfections(
        g, df, interventions, Ns, seedmatrix, i, ngroups, ntimes
    )
end

function samplerenewaldidinfections(g, df, interventions, Ns, seedmatrix, i, ngroups, ntimes)
    _samplerenewaldidinfectionsassertions(df, seedmatrix, i, ngroups, ntimes)
    return _samplerenewaldidinfections(
        g, df, interventions, Ns, seedmatrix, i, ngroups, ntimes
    )
end

function _samplerenewaldidinfections(
    g::Vector, df, interventions, Ns, seedmatrix, i, ngroups, ntimes
)
    return _samplerenewaldidinfections(
        generationtime, df, interventions, Ns, seedmatrix, i, ngroups, ntimes; 
        vec=g,
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


# Assertions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _samplerenewaldidinfectionsassertions(df, seedmatrix, i, ngroups, ntimes)
    _ngroups(seedmatrix) == ngroups || throw(_seedwidthdimensionerror)
    M_xmaxt = _ntimes(seedmatrix) + ntimes
    "M_x[$M_xmaxt, 1]" in names(df) || throw(_dfMxdimensiontoosmallerror())
    "M_x[$(M_xmaxt + 1), 1]" ∉ names(df) || throw(_dfMxdimensiontoolargererror())
    i <= size(df, 1) || throw(BoundsError(df, [i, 1]))
    return nothing
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _dfMxdimensiontoosmallerror()
    return DimensionMismatch(
        "`df` must contain values of M_x for all times up to `_ntimes(seedmatrix) + ntimes`"
    )
end

function _dfMxdimensiontoolargererror()
    return DimensionMismatch(
        "`df` must not contain values M_x for times beyond `_ntimes(seedmatrix) + ntimes`"
    )
end

_seedwidthdimensionerror() = DimensionMismatch("width of `seedmatrix` must equal `ngroups`")
