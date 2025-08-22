# Estimates of generation interval 

# Constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Provided generation intervals

"""
    g_covid(t::Integer)

Return estimate of serial interval for SARS-CoV-2.

`t` is time in days since infection.

Data from https://github.com/mrc-ide/EpiEstim/blob/master/data/covid_deaths_2020_uk.rda

# Examples
```jldoctest
julia> g_covid(10)
0.0382484935

julia> g_covid(0)
0.0
```
"""
function g_covid(t::Integer)
    if t < 0 || t > 30
        return zero(COVIDSERIALINTERVAL[1])
    else
       return COVIDSERIALINTERVAL[t+1]
    end
end

"""
    g_seir(t; mu, kappa=RenewalDiD.automatic)

Estimate generation interval for a susceptible–exposed–infectious–recovered model.

`t` is time in days since infection. `kappa` is the rate of progression from exposed to 
infectious, and `mu` is the recovery rate.

If `mu == kappa` you may input only `mu`.

Equations from D. Champredon, J. Dushoff, and D.J. Earn (2018). Equivalence of the 
Erlang-distributed SEIR epidemic model and the renewal equation. *SIAM J Appl Math*
**78**: 3258–3278.

# Examples
```jldoctest
julia> g_seir(10; mu=0.3)
0.04480836153107755

julia> g_seir(10; mu=0.3) == g_seir(10; mu=0.3, kappa=0.3)
true

julia> g_seir(10; mu=0.3, kappa=0.5)
0.03228684102658386
``` 
"""
function g_seir(t; mu, kappa=automatic)
    if t < 0 
        return zero(_g_seir(zero(t), mu, kappa))
    else
        return _g_seir(t, mu, kappa)
    end
end

_g_seir(t, mu, ::Automatic) = _gseirequalmukappa(t, mu)

function _g_seir(t, mu, kappa)
    mu == kappa && return _gseirequalmukappa(t, mu)
    return _gseirunequalmukappa(t, kappa, mu)
end

_gseirequalmukappa(t, mu) = mu^2 * t * exp(-mu * t)

function _gseirunequalmukappa(t, mu, kappa)
    return kappa * mu * (exp(-mu * t) - exp(-kappa * t)) / (kappa - mu)
end

"""
    vectorg_seir(mu, kappa=automatic; t_max::Integer=28)
    vectorg_seir(; mu, kappa=automatic, t_max::Integer=28)

Produce a vector of generation intervals for a susceptible–exposed–infectious–recovered 
    model.

The length of the vector is `t_max + 1`, i.e. daily from `t=0` to `t=t_max`

`kappa` is the rate of progression from exposed to infectious, and `mu` is the recovery
    rate. These can be entered as positional arguments or keyword arguments.

Equations from D. Champredon, J. Dushoff, and D.J. Earn (2018). Equivalence of the 
Erlang-distributed SEIR epidemic model and the renewal equation. *SIAM J Appl Math*
**78**: 3258–3278.

See also [`g_seir`](@ref).

# Examples
```jldoctest
julia> vectorg_seir(0.4; t_max=5)
6-element Vector{Float64}:
 0.0
 0.10725120736570232
 0.14378526851751092
 0.144573221717857
 0.12921377151657948
 0.10826822658929018

julia> vectorg_seir(0.4, 0.5; t_max=5)
6-element Vector{Float64}:
 0.0
 0.12757877264601186
 0.1628990458915585
 0.15612810352754444
 0.1331224695160854
 0.10650056922542783

julia> vectorg_seir(kappa=0.5, mu=0.4, t_max=5)
6-element Vector{Float64}:
 0.0
 0.12757877264601186
 0.1628990458915585
 0.15612810352754444
 0.1331224695160854
 0.10650056922542783
``` 
"""
vectorg_seir(; mu, kappa=automatic, kwargs...) = vectorg_seir(mu, kappa; kwargs...)

function vectorg_seir(mu, kappa=automatic; t_max::Integer=28)
    v = _vectorg_seir(mu, kappa, t_max)
    for i in eachindex(v)  # remove values of -0.0 for cosmetic reasons 
        if v[i] == 0.0
            v[i] = 0.0 
        end
    end
    return v
end

_vectorg_seir(mu, kappa, t_max) = [_g_seir(t, mu, kappa) for t in 0:1:t_max]

## Functions for user-supplied generation intervals

"""
    generationtime(t::Integer; func=automatic, vec=automatic, <keyword arguments for func>)
    generationtime(func_or_vector, t::Integer; <keyword arguments for func>)
    
Return generation time from a user-supplied vector or function.

If a vector is provided, the first value is taken as `time = 0`.

Keyword arguments are passed to the function provided. 

Note that this function does not check whether the negative values will be returned or if
the sum of outputs will exceed 1. Use `testgenerationtime` for this.

# Examples
```jldoctest
julia> myvec = [0, 0.1, 0.2];

julia> generationtime(0; vec=myvec)
0.0

julia> generationtime(myvec, 1)
0.1

julia> myfunc(t; a) = a * t;

julia> generationtime(1; func=myfunc, a=0.1)
0.1

julia> generationtime(myfunc, 2; a=0.2)
0.4

julia> generationtime(2)
ERROR: ArgumentError: a function or vector must be passed as either a positional or keyword \
    argument
``` 
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
    func == vec && return _generationtime(func, t; kwargs...)
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

"""
    testgenerationtime(func_or_vector; muteinfo, t_max, <keyword arguments for func>)
    
Test suitability of user-supplied vector or function as a generation interval.

Checks that the sum of outputs will be `≤ 1` and that no outputs will be negative. Throws an 
    argument error if either is violated.

# Keyword arguments 
- `muteinfo=false`: whether to provide information about the output if tests passed 
- `t_max=1000`: the maximum time that will be assessed if a function is supplied

Other keyword arguments are passed to the function provided. 

# Examples
```jldoctest
julia> myvec1 = [0, 0.1, 0.2];

julia> testgenerationtime(myvec1)
[ Info: Vector sum is 0.30000000000000004, with minimum value 0.0

julia> testgenerationtime(myvec1; muteinfo=true)

julia> myvec2 = [0, 0.1, 0.2, 0.8];

julia> testgenerationtime(myvec2)
ERROR: ArgumentError: 1.1: sum of all generation times must be ≤ 1

julia> myfunc(t; a) = a * t;

julia> testgenerationtime(myfunc; a=0.1)
ERROR: ArgumentError: 50050.00000000001: sum of all generation times must be ≤ 1 (function \
    tested on x ∈ {0, 1, …, 1000})

julia> testgenerationtime(myfunc; a=0.1, t_max=3)
[ Info: Vector sum is 0.6000000000000001, with minimum value 0.0
``` 
"""
testgenerationtime(x; kwargs...) = _testgenerationtime(x; kwargs...)

function _testgenerationtime(f::Function; t_max::Integer=1000, kwargs...)
    return _testgenerationtime(f, t_max; kwargs...)
end

function _testgenerationtime(f::Function, t_max; muteinfo=false, kwargs...)
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

const _Useablegenerationfunctions = Union{
    typeof(g_covid),
    typeof(g_seir),
    typeof(generationtime)
}


# Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

_g_0_warning(g0) = @warn _g_0_warningtext(g0)

function _g_0_warningtext(g0)
    g0string = "g(0) == $g0"
    remainingstring = ": g(0) is never called and is assumed to equal 0. If you intended \
        this value for g(1) you should recheck the indexing"
    return g0string * remainingstring
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _generationtimebothfunctionvectorerror()
    return ArgumentError("only one of the `func`, `vec` keyword arguments may be used")
end

function _generationtimenofunctionvectorerror()
    m = "a function or vector must be passed as either a positional or keyword argument"
    return ArgumentError(m)
end

function _negativegenerationtimerror(mv, funclengthinfo)
    m = "$mv: generation interval values can never be negative$funclengthinfo"
    return ArgumentError(m)
end

function _toolargegenerationtimerror(sv, funclengthinfo)
    return ArgumentError("$sv: sum of all generation times must be ≤ 1$funclengthinfo")
end
