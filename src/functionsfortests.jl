# functions called by tests 

"""
    RenewalDiD.FittedParameterTestFunctions

Functions that are used when testing the package `RenewalDiD`.

Nothing from this module is exported by the package. The module exports the function 
    `testdataframe`.
"""
module FittedParameterTestFunctions

using DataFrames: DataFrame, insertcols!
using Random: AbstractRNG, default_rng
using RenewalDiD: Automatic, automatic 

export testdataframe

function testdataframe(
    rng::AbstractRNG=default_rng(); 
    nchains, niterations, ngroups, ntimes, nseeds, kwargs...
)
    return testdataframe(rng, nchains, niterations, ngroups, ntimes, nseeds; kwargs...)
end

function testdataframe(nchains, niterations, ngroups, ntimes, nseeds; kwargs...)
    return testdataframe(
        default_rng(), nchains, niterations, ngroups, ntimes, nseeds; 
        kwargs...
    )
end

function testdataframe(
    rng::AbstractRNG, nchains, niterations, ngroups, ntimes, nseeds; 
    gammadefault=automatic, thetadefault=automatic, mxdefault=automatic,
    kwargs...
)
    df = DataFrame(
        :iteration => repeat(1:niterations; outer=nchains),
        :chain => repeat(1:nchains; inner=niterations),
    )
    nrows = nchains * niterations
    _addtotestdataframe!(df, :tau, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :alpha, kwargs, rand(rng, nrows))
    _addtotestdataframe!(df, :sigma_gamma, kwargs, rand(rng, nrows))
    _addgammastotestdataframe!(rng, df, ngroups, nrows, kwargs, gammadefault)
    _addthetastotestdataframe!(rng, df, ntimes, nrows, kwargs, thetadefault)
    _addtotestdataframe!(df, :sigma_theta, kwargs, rand(rng, nrows))
    _addmxtotestdataframe!(rng, df, ngroups, ntimes, nseeds, nrows, kwargs, mxdefault)
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
    insertcols!(
        df,
        name => haskey(kws, name) ? kws[name] : default
    )
    return nothing
end

function _addgammastotestdataframe!(rng, df, ngroups, nrows, kwargs, ::Automatic)
    _addgammastotestdataframe!(rng, df, ngroups, nrows, kwargs, rand(rng, nrows))
    return nothing
end

function _addgammastotestdataframe!(::Any, df, ngroups, ::Any, kwargs, gammadefault)
    for g in 1:(ngroups - 1)
        _addtotestdataframe!(df, Symbol("gammas_raw[$g]"), kwargs, gammadefault)
    end
    return nothing
end

function _addthetastotestdataframe!(rng, df, ntimes, nrows, kwargs, ::Automatic)
    _addthetastotestdataframe!(rng, df, ntimes, nrows, kwargs, rand(rng, nrows))
    return nothing
end

function _addthetastotestdataframe!(::Any, df, ntimes, ::Any, kwargs, thetadefault)
    for t in 1:(ntimes - 1)
        _addtotestdataframe!(df, Symbol("thetas_raw[$t]"), kwargs, thetadefault)
    end
    return nothing
end

function _addmxtotestdataframe!(rng, df, ngroups, ntimes, nseeds, nrows, kwargs, ::Automatic)
    _addmxtotestdataframe!(
        rng, df, ngroups, ntimes, nseeds, nrows, kwargs, rand(rng, nrows)
    )
    return nothing
end

function _addmxtotestdataframe!(::Any, df, ngroups, ntimes, nseeds, ::Any, kwargs, mxdefault)
    for g in 1:(ngroups), t in 1:(ntimes + nseeds)
        _addtotestdataframe!(df, Symbol("M_x[$t, $g]"), kwargs, mxdefault)
    end
    return nothing
end

end  # module FittedParameterTestFunctions
