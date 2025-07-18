# Estimates of generation interval 

## Constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


## Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Provided generation intervals

"""
    g_covid(t::Integer)

Estimate of serial interval for SARS-CoV-2.

`t` is time in days since infection.

Data from https://github.com/mrc-ide/EpiEstim/blob/master/data/covid_deaths_2020_uk.rda
"""
function g_covid(t::Integer)
    if t < 0 || t > 30
        return zero(COVIDSERIALINTERVAL[1])
    else
       return COVIDSERIALINTERVAL[t+1]
    end
end

"""
    g_seir(t; gamma, sigma)

Generation interval for a susceptible–exposed–infectious–recovered model.

`t` is time in days since infection. `sigma` is the rate of progression from exposed to 
infectious, and `gamma` is the recovery rate.

If `gamma == sigma` only input `gamma` and leave `sigma=RenewalDiD.automatic` (the default).

Equations from D. Champredon, J. Dushoff, and D.J. Earn (2018). Equivalence of the 
Erlang-distributed SEIR epidemic model and the renewal equation. *SIAM J Appl Math*
**78**: 3258–3278. 
"""
function g_seir(t; gamma, sigma=automatic)
    if t < 0 
        return zero(g_seir(0; gamma, sigma))
    else
        return _g_seir(t, gamma, sigma)
    end
end

_g_seir(t, gamma, ::Automatic) = _gseirequalgammasigma(t, gamma)

function _g_seir(t, gamma, sigma)
    gamma == sigma && return _gseirequalgammasigma(t, gamma)
    return _gseirunequalgammasigma(t, sigma, gamma)
end

_gseirequalgammasigma(t, gamma) = gamma^2 * t * exp(-gamma * t)

function _gseirunequalgammasigma(t, gamma, sigma)
    return sigma * gamma * (exp(-gamma * t) - exp(-sigma * t)) / (sigma - gamma)
end

function vectorg_seir(gamma, sigma=automatic; t_max::Integer=28)
    return [_g_seir(t, gamma, sigma) for t in 0:1:t_max]
end

## Functions for user-supplied generation intervals

"""
    generationtime(t::Integer; func=automatic, vec=automatic, kwargs...)
    generationtime(f_or_v, t::Integer; kwargs...)
    
Return generation time from a user-supplied vector or function.

Note that this function does not check whether the negative values will be returned or if
the sum of outputs will exceed 1. Use `testgenerationtime` for this.
"""
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
    throw(_generationtimenofunctionvectorerror())
    return nothing
end

function _generationtime(func, vec, ::Integer; kwargs...)
    throw(_generationtimebothfunctionvectorerror())
    return nothing
end

function _generationtime(f::Function, t::Integer; t_max=automatic, kwargs...)
    return __generationtime(f, t, t_max; kwargs...)
end

function _generationtime(v::AbstractVector, t::Integer)
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

testgenerationtime(x; kwargs...) = _testgenerationtime(x; kwargs...)

function _testgenerationtime(f::Function; t_max=automatic, kwargs...)
    return _testgenerationtime(f, t_max; kwargs...)
end

function _testgenerationtime(f::Function, ::Automatic; kwargs...)
    return _testgenerationtime(f, 1000; kwargs...)
end

function _testgenerationtime(f::Function, t_max::Number; muteinfo=false, kwargs...)
    v = [f(x; kwargs...) for x in 0:1:t_max]
    return _testgenerationtime(
        v; 
        funclengthinfo=" (function tested on x ∈ {0, 1, …, $t_max})", muteinfo
    )
end

function _testgenerationtime(v::AbstractVector; funclengthinfo="", muteinfo=false)
    sv = sum(v)
    sv <= 1  || throw(_toolargegenerationtimerror(sv, funclengthinfo))
    mv = minimum(v) 
    mv >= 0 || throw(_negativegenerationtimerror(mv, funclengthinfo))
    sum(v) <= 1
    v[1] == 0 || _g_0_warning(v[1])
    muteinfo || @info "Vector sum is $sv, with minimum value $mv"
    return nothing 
end 


# Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

_g_0_warning(g0) = @warn _g_0_warningtext(g0)

function _g_0_warningtext(g0)
    g0string = "g(0) == $g0"
    remainingstring = ": g(0) is never called and is assumed to equal 0. If you intended \
        this value for g(1) you should recheck the indexing"
    return g0string * remainingstring
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _g_seirequalgammasigmaerror(gamma, sigma)
    return ArgumentError(
        "$gamma, $sigma: if gamma==sigma only the gamma keyword argument should be provided"
    )
end

function _generationtimebothfunctionvectorerror()
    return ArgumentError("only one of the `func`, `vec` keyword arguments may be used")
end

function _generationtimenofunctionvectorerror()
    return ArgumentError("""
        a function or vector must be passed as either a positional or keyword argument
    """)
end

function _negativegenerationtimerror(mv, funclengthinfo)
    return ArgumentError("$mv: generation interval values can never be negative$funclengthinfo")
end

function _toolargegenerationtimerror(sv, funclengthinfo)
    return ArgumentError("$sv: sum of all generation times must be ≤ 1$funclengthinfo")
end
