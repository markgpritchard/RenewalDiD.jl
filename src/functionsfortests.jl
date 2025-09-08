# functions that are used when testing the package `RenewalDiD`

# nothing from this file is exported

function testdataframe(
    rng::AbstractRNG=default_rng(); 
    nchains, niterations, ngroups, ntimes, ninterventions=1, nseeds, kwargs...
)  
    # default `ninterventions` given as this was not previously an option
    return testdataframe(
        rng, nchains, niterations, ngroups, ntimes, ninterventions, nseeds; 
        kwargs...
    )
end

function testdataframe(
    nchains, niterations, ngroups, ntimes, ninterventions, nseeds; 
    kwargs...
)
    return testdataframe(
        default_rng(), nchains, niterations, ngroups, ntimes, ninterventions, nseeds; 
        kwargs...
    )
end

function testdataframe(
    rng::AbstractRNG, nchains, niterations, ngroups, ntimes, ninterventions, nseeds; 
    gammadefault=automatic, 
    thetadefault=automatic, 
    taudefault=automatic,
    mxdefault=automatic, 
    predictiondefault=automatic,
    kwargs...
)
    df = DataFrame(
        :iteration => repeat(1:niterations; outer=nchains),
        :chain => repeat(1:nchains; inner=niterations),
    )
    nrows = nchains * niterations
    _addtaustotestdataframe!(rng, df, ninterventions, nrows, kwargs, taudefault)
    _addtotestdataframe!(df, :alpha, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :sigma_gamma, kwargs, rand(rng, nrows))
    _addgammastotestdataframe!(rng, df, ngroups, nrows, kwargs, gammadefault)
    _addthetastotestdataframe!(rng, df, ntimes, nrows, kwargs, thetadefault)
    _addtotestdataframe!(df, :sigma_theta, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :psi, kwargs, rand(rng, nrows))
    _addmxtotestdataframe!(rng, df, ngroups, ntimes, nseeds, nrows, kwargs, mxdefault)
    _addtotestdataframe!(df, :minsigma2, kwargs, rand(rng, Beta(1, 2), nrows))
    _addtotestdataframe!(df, :fittingsigma, kwargs, -100 .* rand(rng, nrows))
    _addpredobsinfsigmatotestdataframe!(
        rng, df, ngroups, ntimes, nrows, kwargs, predictiondefault
    )
    _addtotestdataframe!(df, :lp, kwargs, -100 .* rand(rng, nrows))
    _addtotestdataframe!(df, :n_steps, kwargs, ones(nrows))
    _addtotestdataframe!(df, :is_accept, kwargs, ones(nrows))
    _addtotestdataframe!(df, :acceptance_rate, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :log_density, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :hamiltonian_energy, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :hamiltonian_energy_error, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :max_hamiltonian_energy_error, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :tree_depth, kwargs, 3 .* ones(nrows))
    _addtotestdataframe!(df, :numerical_error, kwargs, zeros(nrows))
    _addtotestdataframe!(df, :step_size, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :nom_step_size, kwargs, rand(rng, nrows))
    return df
end

function _addtotestdataframe!(df, name, kws, default)
    return insertcols!(
        df,
        name => haskey(kws, name) ? kws[name] : default
    )
end

function _addgammastotestdataframe!(rng::AbstractRNG, df, ngroups, nrows, kws, ::Automatic)
    return _addgammastotestdataframe!(rng, df, ngroups, nrows, kws, rand(rng, nrows))
end

function _addgammastotestdataframe!(
    ::AbstractRNG, df, ngroups, ::Int, kws, gammavalues::AbstractVector
)
    for g in 1:(ngroups - 1)
        _addtotestdataframe!(df, Symbol("gammas_raw[$g]"), kws, gammavalues)
    end
    return nothing
end

function _addthetastotestdataframe!(rng::AbstractRNG, df, ntimes, nrows, kws, ::Automatic)
    return _addthetastotestdataframe!(rng, df, ntimes, nrows, kws, rand(rng, nrows))
end

function _addthetastotestdataframe!(
    ::AbstractRNG, df, ntimes, ::Int, kws, thetavalues::AbstractVector
)
    for t in 1:(ntimes - 1)
        _addtotestdataframe!(df, Symbol("thetas_raw[$t]"), kws, thetavalues)
    end
    return nothing
end

function _addtaustotestdataframe!(
    rng::AbstractRNG, df, ninterventions, nrows, kws, ::Automatic
)
    return _addtaustotestdataframe!(rng, df, ninterventions, nrows, kws, rand(rng, nrows))
end

function _addtaustotestdataframe!(
    ::AbstractRNG, df, ninterventions, ::Int, kws, tauvalues::AbstractVector
)
    for k in 1:ninterventions
        _addtotestdataframe!(df, Symbol("tau[$k]"), kws, tauvalues)
    end
    return nothing
end

function _addmxtotestdataframe!(
    rng::AbstractRNG, df, ngroups, ntimes, nseeds, nrows, kws, ::Automatic
)
    return _addmxtotestdataframe!(
        rng, df, ngroups, ntimes, nseeds, nrows, kws, rand(rng, nrows)
    )
end

function _addmxtotestdataframe!(
    ::AbstractRNG, df, ngroups, ntimes, nseeds, ::Int, kws, mxvalues::AbstractVector
)
    for g in 1:(ngroups), t in 1:(ntimes + nseeds)
        _addtotestdataframe!(df, Symbol("M_x[$t, $g]"), kws, mxvalues)
    end
    return nothing
end

function _addpredobsinfsigmatotestdataframe!(
    rng::AbstractRNG, df, ngroups, ntimes, nrows, kws, ::Automatic
)
    return _addpredobsinfsigmatotestdataframe!(
        rng, df, ngroups, ntimes, nrows, kws, rand(rng, nrows)
    )
end

function _addpredobsinfsigmatotestdataframe!(
    ::AbstractRNG, df, ngroups, ntimes, ::Int, kws, predictionvalues::AbstractVector
)
    for g in 1:(ngroups), t in 1:(ntimes + 1)
        _addtotestdataframe!(
            df, 
            Symbol("predictobservedinfectionssigmamatrix[$t, $g]"), 
            kws, 
            predictionvalues
        )
    end
    return nothing
end

function testsimulation(rng::AbstractRNG=default_rng())
    u0_1 = simulationu0(; s=98, e=2)
    u0_2 = simulationu0(; s=198, e=2)
    u0_3 = simulationu0(; s=48, e=2)
    mu = 0.2
    delta = 0.3
    psi = 0.6
    kappa = 0.5
    s1 = packsimulationtuple( ; 
        u0=u0_1, beta=_testsimulationbeta1, mu, delta, psi, kappa, intervention=nothing,
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=_testsimulationbeta2, mu, delta, psi, kappa, intervention=4,
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=_testsimulationbeta3, mu, delta, psi, kappa, intervention=6,
    )
    return packsimulations(rng, 10, s1, s2, s3)
end

_testsimulationbeta1(t) = 0.3 + 0.1 * cos((t - 20) * 2pi / 365)
_testsimulationbeta2(t) = _testsimulationbeta1(t) * t < 4 ? 1.1 : 0.88
_testsimulationbeta3(t) = _testsimulationbeta1(t) * t < 6 ? 0.9 : 0.72

function tupleforsamplerenewaldidinfections(data; vec=nothing, kwargs...)
    return _tupleforsamplerenewaldidinfections(data, vec; kwargs...)
end

function _tupleforsamplerenewaldidinfections(d, vec::AbstractVector; delaydistn=Normal(0, 0))
    return Model(
        _renewaldid,
        (
            observedcases=_observedcases(d),
            interventions=_interventions(d),
            expectedseedcases=_expectedseedcases(d),
            Ns=_ns(d),
            g=generationtime,
            delaydistn=delaydistn,
            n_seeds=_nseeds(d),
        ),
        (vec=vec, ),
    )
end
