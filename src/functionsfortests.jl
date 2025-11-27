# functions that are used when testing the package `RenewalDiD`

"""
    RenewalDiD.testdataframe([rng,] nchains, niterations, ngroups, ntimes, ninterventions, nseeds; <keyword arguments)
    RenewalDiD.testdataframe([rng]; nchains, niterations, ngroups, ntimes, ninterventions, nseeds, <keyword arguments)

Generate an example `DataFrame` used in tests.

# Arguments 
- `rng::AbstractRNG=default_rng()`: random number generator
- `nchains`: number of chains 
- `niterations`: number of iterations per chain
- `ngroups`: number of groups
- `ntimes`: number of times
- `ninterventions=1`: number of interventions 
- `nseeds`: number of rows used to seed the renewal equation
Other keyword arguments supply default values, either for individual columns in the DataFrame or groups:
- `gammadefault=RenewalDiD.automatic` 
- `thetadefault=RenewalDiD.automatic` 
- `taudefault=RenewalDiD.automatic`
- `mxdefault=RenewalDiD.automatic`
- `predictiondefault=RenewalDiD.automatic`

It is expected that an odd number of quantiles will be supplied, symmetrical around the 
    median (`0.5`)

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(100);

julia> RenewalDiD.testdataframe(
       rng;
       nchains=3, niterations=1, ngroups=1, ntimes=1, nseeds=3, alpha=1, psi=0.2,      
       )
3×27 DataFrame
 Row │ iteration  chain  tau[1]     alpha  sigma_gamma  sigma_theta  psi       ⋯
     │ Int64      Int64  Float64    Int64  Float64      Float64      Float64   ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │         1      1  0.208041       1     0.548821     0.22052       0.2   ⋯
   2 │         1      2  0.181642       1     0.371745     0.371992      0.2
   3 │         1      3  0.0140921      1     0.402931     0.262554      0.2
                                                              20 columns omitted
```
"""
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
        _addtotestdataframe!(df, Symbol("logtau[$k]"), kws, tauvalues)
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

"""
    RenewalDiD.testsimulation([rng])

Generate a simulation.

This function uses a default set of arguments and is intended for testing other functions 
    rather than as a simulation to analyse.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(10);

julia> RenewalDiD.testsimulation(rng)
SimulationData{Int64, InterventionMatrix{Int64}, Vector{Int64}}
 observedcases:  [0 0 0; 0 0 0; … ; 0 2 2; 1 7 1]
 interventions:  [0 0 0; 0 0 0; … ; 0 1 1; 0 1 1] {duration 10, starttimes [nothing, 4, 6]}
 Ns:             [100, 200, 50]
 exptdseedcases: [0.0 0.0 0.0; 0.0 0.0 0.0; … ; 0.0 0.0 0.0; 0.5 0.5 0.5]
```
"""
function testsimulation(rng::AbstractRNG=default_rng())
    u0_1 = simulationu0(; S=98, E=2)
    u0_2 = simulationu0(; S=198, E=2)
    u0_3 = simulationu0(; S=48, E=2)
    eta = 0.2
    phi = 0.6
    sigma = 0.5
    s1 = packsimulationtuple( ; 
        u0=u0_1, beta=_testsimulationbeta1, sigma, eta, phi, intervention=nothing,
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=_testsimulationbeta2, sigma, eta, phi, intervention=4,
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=_testsimulationbeta3, sigma, eta, phi, intervention=6,
    )
    return packsimulations(rng, 10, s1, s2, s3)
end

_testsimulationbeta1(t) = 0.3 + 0.1 * cos((t - 20) * 2pi / 365)
_testsimulationbeta2(t) = _testsimulationbeta1(t) * t < 4 ? 1.1 : 0.88
_testsimulationbeta3(t) = _testsimulationbeta1(t) * t < 6 ? 0.9 : 0.72

"""
    RenewalDiD.testmodel(data::AbstractRenewalDiDData; vec::AbstractVector)

Generate a `renewaldid` model for use in tests.
    
Uses `generationtime` for the generation interval, with a vector `vec`.

# Examples
```jldoctest
julia> data = RenewalDiDData( ;
       observedcases=zeros(11, 3),
       interventions=zeros(10, 3),
       exptdseedcases=zeros(7, 3),
       );

julia> RenewalDiD.testmodel(data; vec=zeros(2))
DynamicPPL.Model{typeof(RenewalDiD._renewaldid), (:observedcases, :interventions, \
    :expectedseedcases, :Ns, :g, :delaydistn, :n_seeds), (:vec,), (), \
    Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Nothing, \
    typeof(generationtime), Normal{Float64}, Int64}, Tuple{Vector{Float64}}, \
    DynamicPPL.DefaultContext}(RenewalDiD._renewaldid, (observedcases = [0.0 0.0 0.0; 0.0 \
    0.0 0.0; … ; 0.0 0.0 0.0; 0.0 0.0 0.0], interventions = [0.0 0.0 0.0; 0.0 0.0 0.0; … ; \
    0.0 0.0 0.0; 0.0 0.0 0.0], expectedseedcases = [0.0 0.0 0.0; 0.0 0.0 0.0; … ; 0.0 0.0 \
    0.0; 0.0 0.0 0.0], Ns = nothing, g = generationtime, delaydistn = \
    Normal{Float64}(μ=0.0, σ=0.0), n_seeds = 7), (vec = [0.0, 0.0],), \
    DynamicPPL.DefaultContext())
```
"""
function testmodel(data; vec=nothing, kwargs...)
    return _testmodel(data, vec; kwargs...)
end

function _testmodel(d, vec::AbstractVector; delaydistn=Normal(0, 0))
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
