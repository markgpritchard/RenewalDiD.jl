const COVIDSERIALINTERVAL = [  # not exported but can be accessed via `g_covid`
# result from https://github.com/mrc-ide/EpiEstim/blob/master/data/covid_deaths_2020_uk.rda
    0.0000000000,
    0.0440204506,
    0.1298284450,
    0.1397552873,
    0.1277718301,
    0.1100166556,
    0.0917470443,
    0.0749977679,
    0.0604725660,
    0.0482765015,
    0.0382484935,
    0.0301228893,
    0.0236092441,
    0.0184305583,
    0.0143398489,
    0.0111254255,
    0.0086104507,
    0.0066498251,
    0.0051260438,
    0.0039448946,
    0.0030314300,
    0.0023264019,
    0.0017832132,
    0.0013653739,
    0.0010444113,
    0.0007981781,
    0.0006094926,
    0.0004650564,
    0.0003545982,
    0.0002701988,
    0.0002057625
]

function g_covid(t::Integer)
    if t < 0 || t > 30
        return zero(COVIDSERIALINTERVAL[1])
    else
       return COVIDSERIALINTERVAL[t+1]
    end
end

g_seir(t; gamma, sigma=automatic) = _g_seir(t, gamma, sigma)

function _g_seir(t, gamma, sigma)
    if t < 0 
        return __g_seir(0, gamma, sigma)
    else
        return __g_seir(t, gamma, sigma) 
    end
end

function __g_seir(t, gamma, sigma)
    if gamma == sigma 
        return _gseirequalgammasigma(t, gamma)
    else
        return _gseirunequalgammasigma(t, sigma, gamma)
    end
end

__g_seir(t, gamma, ::Automatic) = _gseirequalgammasigma(t, gamma)

_gseirequalgammasigma(t, gamma) = gamma^2 * t * exp(-gamma * t)

function _gseirunequalgammasigma(t, gamma, sigma)
    return sigma * gamma * (exp(-gamma * t) - exp(-sigma * t)) / (sigma - gamma)
end

function vectorg_seir(gamma, sigma=automatic; t_max::Integer=28)
    return [_g_seir(t, gamma, sigma) for t in 0:1:t_max]
end

function generationtime(t::Integer; func=automatic, vec=automatic, kwargs...)
    return _generationtime(func, vec, t; kwargs...)
end

generationtime(f_or_v, t::Integer; kwargs...) = _generationtime(f_or_v, t; kwargs...)

function _generationtime(func, ::Automatic, t::Integer; kwargs...)
    return _generationtime(func, t; kwargs...)
end

function _generationtime(::Automatic, vec, t::Integer; kwargs...)
    return _generationtime(vec, t; kwargs...)
end

function _generationtime(::Automatic, ::Automatic, ::Integer; kwargs...)
    throw(ArgumentError("""
        a function or vector must be passed as either a positional or keyword argument
    """))
end

function _generationtime(func, vec, ::Integer; kwargs...)
    throw(ArgumentError("only one of the `func`, `vec` keyword arguments may be used"))
end

function _generationtime(f::Function, t::Integer; t_max=automatic, kwargs...)
    _testgenerationtime(f, t_max; kwargs...)
    return __generationtime(f, t, t_max; kwargs...)
end

function _generationtime(v::AbstractVector, t::Integer)
    _testgenerationtime(v)
    if t < 0 || t >= length(v)
        return zero(v[1])
    else
        return v[t+1]
    end
end

function __generationtime(f::Function, t::Integer, ::Automatic; kwargs...)
    if t < 0 
        return zero(f(0; kwargs...))
    else
        return f(t; kwargs...)
    end
end

function __generationtime(f::Function, t::Integer, t_max::Integer; kwargs...)
    if t < 0 || t > t_max
        return zero(f(0; kwargs...))
    else
        return f(t; kwargs...)
    end
end

function _testgenerationtime(v::AbstractVector)
    return __testgenerationtime(
        v,
        "all values of v must be positive",
        "sum of v must be ≤ 1"
    )
end

function _testgenerationtime(f::Function, ::Automatic; kwargs...)
    return _testgenerationtime(f, 1000; kwargs...)
end

function _testgenerationtime(f::Function, t_max::Integer; kwargs...)
    v = [f(x; kwargs...) for x in 0:1:t_max]
    return __testgenerationtime(
        v,
        "all values of f(x) must be positive (tested on x ∈ {0, 1, …, $t_max})",
        "sum of f(x) must be ≤ 1 (tested on x ∈ {0, 1, …, $t_max})"
    )
end

function __testgenerationtime(v, minerror, sumerror)
    minimum(v) >= 0 || throw(ArgumentError("$(minimum(v)): $minerror"))
    sum(v) <= 1 || throw(ArgumentError("$(sum(v)): $sumerror"))
    return nothing
end
