# generates simulated datasets

# `runsimulation` is the main function in this file. All parameters and compartments have
#  the meanings described in this function's documentation.
"""
    runsimulation([rng], duration, u0; beta, gamma, delta, theta, sigma)
    runsimulation([rng], duration, u0, beta, gamma, delta, theta, sigma)

Generate simulated data.

Model uses Gillespie's stochastic continuous time method with 6 compartments:
- `s`: susceptible 
- `e`: exposed (infected but not yet infectious)
- `i_n`: infectious and will never be diagnosed 
- `i_f`: infectious and will be diagnosed before recovery 
- `i_d`: infectious and diagnosed 
- `r`: recovered

A 7th compartment in the model records the cumulative number of diagnoses.

Parameters `beta`, `gamma`, `delta` `theta` and `sigma` may be numbers, or may be functions 
taking time as their only argument to give time-varying values. They can be entered as 
positional arguments or keyword arguments. All parameters must be non-negative.

# Arguments
- `rng::AbstractRNG=Random.default_rng()`: random number generator
- `duration::Integer`: duration of simulation
- `u0::AbstractVector{<:Integer}`: initial conditions, must be of length 7 describing
    compartments in the order above
- `beta`: transmission parameter
- `gamma`: recovery rate
- `delta`: diagnosis rate (i.e. for those who will be diagnosed, how quickly does it 
    happen?); must satisfy `delta > gamma` (i.e. must be diagnosed before the end of the 
    infectious period)
- `theta`: proportion diagnosed; must satisfy `0 ≤ theta ≤ 1`
- `sigma`: rate of progression from exposed to infectious

# Returns
Returns a matrix of height `duration + 1` giving numbers in each compartment at the end of 
each day. The first row is the conditions supplied in `u0`.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(1);

julia> u0 = simulationu0(; s=100, e=5, n=200);

julia> runsimulation(rng, 10, u0; beta=0.6, gamma=0.25, delta=0.33, theta=0.5, sigma=0.4)
11×7 Matrix{Int64}:
 100  5  0  0  0  95  0
 100  5  0  0  0  95  0
 100  4  0  1  0  95  0
  99  4  0  2  0  95  0
  99  3  0  3  0  95  0
  99  2  0  3  0  96  1
  98  1  1  3  0  97  1
  97  1  1  4  0  97  1
  95  2  2  3  1  97  2
  93  2  4  3  1  97  2
  92  0  5  4  0  99  3

julia> mybetafunc(t) = t <= 5 ? 1 : 0.1;

julia> runsimulation(rng, 10, u0;
       beta=mybetafunc, gamma=0.25, delta=0.33, theta=0.5, sigma=0.4)
11×7 Matrix{Int64}:
 100  5  0  0  0   95  0
  98  5  0  2  0   95  0
  97  3  1  1  1   97  2
  96  3  1  2  1   97  2
  92  5  1  2  3   97  4
  92  3  1  3  1  100  5
  92  2  2  3  0  101  5
  92  2  1  3  0  102  5
  92  1  1  3  1  102  6
  92  0  1  4  1  102  6
  92  0  1  3  1  103  7
```
"""
function runsimulation end

# Constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const _SEIREVENTSMATRIX = [  # matrix of movements between compartments 
    # s    e    i_n  i_f  i_d  r    cumulativediagnoses
     -1    1    0    0    0    0    0  # infections
      0   -1    1    0    0    0    0  # disease progression, subset who won't be diagnosed
      0   -1    0    1    0    0    0  # disease progression, subset who will be diagnosed
      0    0    0   -1    1    0    1  # diagnosis
      0    0   -1    0    0    1    0  # recovery from i_n
      0    0    0    0   -1    1    0  # recovery from i_d
]


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Simulation u0

"""
    simulationu0(; s, e, i_n, i_d, i_f, r, n)

Produce vector of initial conditions for simulation.

# Keyword arguments 
The following represent the model compartments, 
- `s::Integer=0`: susceptible 
- `e::Integer=0`: exposed (infected but not yet infectious)
- `i_n::Integer=0`: infectious and will never be diagnosed 
- `i_f::Integer=0`: infectious and will be diagnosed before recovery 
- `i_d::Integer=0`: infectious and diagnosed 
- `r::Union{<:Integer, Automatic}=automatic`: recovered

The argument `n::Union{<:Integer, Nothing}=nothing` is the population size. If `n` is 
    provided and `r` is not then the remaining population not assigned to other compartments 
    will be assigned to `r`. An error is thrown if a unique non-negative value of `r` cannot
    be calculated.

All arguments are optional.

# Examples
```jldoctest
julia> simulationu0(; s=99, e=1)
7-element Vector{Int64}:
 99
  1
  0
  0
  0
  0
  0

julia> simulationu0(; s=99, e=1, n=200)
7-element Vector{Int64}:
  99
   1
   0
   0
   0
 100
   0

julia> simulationu0(; s=99, e=1, r=100, n=250)
ERROR: ArgumentError: Inconsistent values of `n` and `r` supplied. Calculated n=200 but \
keyword argument n=250.
```
"""
function simulationu0( ; 
    s::Integer=0, 
    e::Integer=0, 
    i_n::Integer=0, 
    i_d::Integer=0, 
    i_f::Integer=0, 
    r::Union{<:Integer, Automatic}=automatic,
    n::Union{<:Integer, Nothing}=nothing,
)
    return _simulationu0(s, e, i_n, i_f, i_d, r, n)
end

function _simulationu0(s, e, i_n, i_f, i_d, ::Automatic, ::Nothing)
    return _simulationu0(s, e, i_n, i_f, i_d, 0, nothing)
end

function _simulationu0(s, e, i_n, i_f, i_d, ::Automatic, n::Integer)
    r = n - (s + e + i_n + i_f + i_d)
    r >= 0 || throw(_simulationu0_n_error(n, n - r))
    return _simulationu0(s, e, i_n, i_f, i_d, r, nothing)
end

function _simulationu0(s, e, i_n, i_f, i_d, r::Integer, n::Integer)
    s + e + i_n + i_f + i_d + r == n || throw(
        _simulationu0_nr_error(n, s + e + i_n + i_f + i_d + r)
    )
    return _simulationu0(s, e, i_n, i_f, i_d, r, nothing)
end

function _simulationu0(s, e, i_n, i_f, i_d, r::Integer, ::Nothing)
    u0 = [s, e, i_n, i_f, i_d, r, 0]
    minimum(u0) >= 0 || throw(_simulationu0_negativeerror(minimum(u0)))
    return u0
end

## Parameters

# Model parameters may be functions or numbers. The `_parameter` function allows either to 
# be handled. In all cases the value must be >= 0. `symbol` is shown error messages when 
# outputs are not as expected. 

# generic parameter name: prefer to use a symbol that is the name of a specific parameter
_parameter(x, t; kwargs...) = _parameter(x, t, :parameter; kwargs...)

_parameter(f::Function, t, symbol; kwargs...) = _parameter(f(t), t, symbol; kwargs...)

function _parameter(x::Number, t, symbol; upper=nothing)
    x >= 0 || throw(_parameterzeroerrorexception(x, t, symbol))
    return __parameter(x, t, symbol, upper)
end

__parameter(x, ::Any, ::Any, ::Nothing) = x

function __parameter(x, t, symbol, upper)
    x <= upper || throw(_parametermaximumerrorexception(x, t, symbol, upper))
    return x
end

_pbeta(beta, t) = _parameter(beta, t, :beta)
_pdelta(delta, t) = _parameter(delta, t, :delta)
_pgamma(gamma, t) = _parameter(gamma, t, :gamma)
_psigma(sigma, t) = _parameter(sigma, t, :sigma)
_ptheta(theta, t) = _parameter(theta, t, :theta; upper=1)  # theta is a proportion 0 ≤ θ ≤ 1

## Force of infection 

# i_n, i_f and i_d assumed to be equally infectious
_foi(beta, i_n, i_f, i_d, n, t) = _pbeta(beta, t) * (i_n + i_f + i_d) / n

## Event rates

_simulatedinfections(beta, s, i_n, i_f, i_d, n, t) = s * _foi(beta, i_n, i_f, i_d, n, t)

function _diseaseprogression_to_i_n(x, theta, sigma, t)
    return x * (1 - _ptheta(theta, t)) * _psigma(sigma, t)
end

_diseaseprogression_to_i_f(x, theta, sigma, t) = x * _ptheta(theta, t) * _psigma(sigma, t)
_diagnosis(x, delta, t) = x * _pdelta(delta, t)
_recovery(x, gamma, t) = x * _pgamma(gamma, t)

# diagnosed individuals have already been infectious for a period represented by `1/delta`
# so recover after a period `1/gamma - 1/delta`
function _recovery(x, gamma, delta, t) 
    timetodiagnosis = 1 / _pdelta(delta, t)
    timetorecovery = 1 / _pgamma(gamma, t)
    if isinf(timetorecovery)
        # no need for error message about time taken to diagnosis if recovery never happens
        remainingtime = timetorecovery
    else 
        timetodiagnosis < timetorecovery || throw(
            _recoveryerrorexception(timetodiagnosis, timetorecovery)
        )
        remainingtime = timetorecovery - timetodiagnosis
    end
    remainingrate = 1 / remainingtime
    return x * remainingrate
end

function _n_seir(u::AbstractVector{<:Integer})
    length(u) == 7 || throw(_ulengtherror(u))
    minimum(u) >= 0 || throw(_negativeuerror(u))
    return sum(@view u[1:6])  # s, e, i_n, i_f, i_d, r 
end

function _seirrates(u::AbstractVector{<:Integer}, t, beta, gamma, delta, theta, sigma)
    n = _n_seir(u)
    s, e, i_n, i_f, i_d, = u
    return [
        _simulatedinfections(beta, s, i_n, i_f, i_d, n, t),  # infections
        _diseaseprogression_to_i_n(e, theta, sigma, t),  # disease progression, not diagnosed
        _diseaseprogression_to_i_f(e, theta, sigma, t),  # disease progression, to be diagnosed
        _diagnosis(i_f, delta, t),  # diagnosis
        _recovery(i_n, gamma, t),  # recovery from i_n
        _recovery(i_d, gamma, delta, t)  # recovery from i_d
    ]
end

## Effect the next event

_tstep(rng, rates) = -log(rand(rng)) / sum(rates)
_nextevent(rng, rates) = sample(rng, eachindex(rates), Weights(rates))
_updateevent!(u, nextevent) = u .+= _SEIREVENTSMATRIX[nextevent, :]

## Simulate a day

_simulateday!(args...) = _simulateday!(default_rng(), args...)

function _simulateday!(rng::AbstractRNG, u, t, beta, gamma, delta, theta, sigma)
    nextday = t + 1 
    while t < nextday
        rates = _seirrates(u, t, beta, gamma, delta, theta, sigma)
        nextevent = _nextevent(rng, rates)
        t += _tstep(rng, rates)
        if t < nextday 
            _updateevent!(u, nextevent) 
        end
    end
    return u 
end

## Whole simulation

runsimulation(args...; kwargs...) = runsimulation(default_rng(), args...; kwargs...) 

function runsimulation(rng::AbstractRNG, duration, u0; beta, gamma, delta, theta, sigma)
    return _runsimulation(rng, duration, u0, beta, gamma, delta, theta, sigma)
end

function runsimulation(rng::AbstractRNG, duration, u0, beta, gamma, delta, theta, sigma)
    return _runsimulation(rng, duration, u0, beta, gamma, delta, theta, sigma) 
end

function _runsimulation(rng, duration, u0, beta, gamma, delta, theta, sigma)
    duration >= 1 || throw(_negativedurationerror(duration))
    output = zeros(Int, duration + 1, 7)
    u = deepcopy(u0)
    output[1, :] = u
    _runsimulationdays!(rng, output, duration, u, beta, gamma, delta, theta, sigma)
    return output 
end

function _runsimulationdays!(rng, output, duration, u, beta, gamma, delta, theta, sigma)
    for t in 1:duration 
        _simulateday!(rng, u, t, beta, gamma, delta, theta, sigma)
        output[t+1, :] = u
    end
    return nothing
end

"""
    simulationcases(M::Matrix)
    simulationcases([rng], duration, u0, beta, gamma, delta, theta, sigma)
    simulationcases([rng], duration, u0; beta, gamma, delta, theta, sigma)

Produce a vector of simulated diagnosed infections.

Either takes a matrix that is the output of `runsimulation` or a set of arguments that gets 
passed to `runsimulation`. 

See `runsimulation` for more details.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng1 = StableRNG(1);

julia> rng2 = StableRNG(1);

julia> u0 = simulationu0(; s=100, e=5, n=200);

julia> runsimulation(rng1, 10, u0;
       beta=0.6, gamma=0.25, delta=0.33, theta=0.5, sigma=0.4);

julia> simulationcases(sim)
11-element Vector{Int64}:
 0
 0
 0
 0
 0
 1
 0
 0
 1
 0
 1

julia> simulationcases(rng2, 10, u0; beta=0.6, gamma=0.25, delta=0.33, theta=0.5, sigma=0.4)
11-element Vector{Int64}:
 0
 0
 0
 0
 0
 1
 0
 0
 1
 0
 1
```
"""
function simulationcases(M::Matrix)
    size(M, 2) == 7 || throw(_simmatrixwidtherror(M))
    duration = size(M, 1) - 1
    return _simulationcases(M, duration)
end

simulationcases(args...; kwargs...) = simulationcases(default_rng(), args...; kwargs...)

function simulationcases(rng::AbstractRNG, duration, args...; kwargs...)
    M = runsimulation(rng, duration, args...; kwargs...)
    return _simulationcases(M, duration)
end

function _simulationcases(M::Matrix{T}, duration) where T
    xs = zeros(T, duration + 1)
    xs[1] == M[1, 7]
    for t in 2:(duration + 1)
        xs[t] = M[t, 7] - M[t-1, 7]
    end
    return xs
end

## Group / pack simulations 

"""
    packsimulations([rng::AbstractRNG], duration, m1_args, args...; <keyword arguments>)

Run a series of simulations and collate results into a `RenewalDiDData` struct to be passed 
    to `renewaldid`.

# Arguments

- `duration`: duration of the simulation.
- `m1_args`: a `Tuple` containing `{u0, beta, gamma, delta, theta, sigma, intervention}` in 
    order for the first group's simulation. The function `packsimulationtuple` can be used 
    to generate this with keyword arguments
- `args...`: equivalent Tuples for remaining groups.
- keyword arguments are passed to `RenewalDiDData`

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(100);

julia> u0_1 = simulationu0(; s=100_000, e=5, i_n=3, i_f=2);

julia> u0_2 = simulationu0(; s=200_000, e=5, i_n=3, i_f=2);

julia> gamma = 0.2; delta = 0.3; theta = 0.6; sigma = 0.5;

julia> _beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365);

julia> _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t);

julia> _beta2(t) = 0.72 * _beta1counter(t);

julia> s1 = packsimulationtuple( ;
       u0=u0_1, beta=_beta1, gamma, delta, theta, sigma, intervention=50,
       );

julia> s2 = packsimulationtuple( ;
       u0=u0_2, beta=_beta2, gamma, delta, theta, sigma, intervention=nothing,
       );

julia> packsimulations(rng, 100, s1, s2)
RenewalDiDData{Int64, InterventionMatrix{Int64}}
 observedcases:  [0 0; 0 0; … ; 6 1380; 9 1414]
 interventions:  [0 0; 0 0; … ; 1 0; 1 0] {duration 100, starttimes [50, nothing]}
 Ns:             [100010, 200010]
 exptdseedcases: [0.0 0.46548963316975533; 0.024913053640556182 0.0; … ; \
    0.4209971568176677 0.0; 0.5200181826119457 0.03451036683024468]
```
"""
function packsimulations(args...; kwargs...)
    return _packsimulations(RenewalDiDData, args...; kwargs...)
end

"""
    packsimulationsunlimitedpopulation([rng::AbstractRNG], duration, m1_args, args...; \
        <keyword arguments>)

Run a series of simulations and collate results into a `RenewalDiDDataUnlimitedPopn` struct 
    to be passed to `renewaldid`.

# Arguments

- `duration`: duration of the simulation.
- `m1_args`: a `Tuple` containing `{u0, beta, gamma, delta, theta, sigma, intervention}` in 
    order for the first group's simulation. The function `packsimulationtuple` can be used 
    to generate this with keyword arguments
- `args...`: equivalent Tuples for remaining groups.
- keyword arguments are passed to `RenewalDiDDataUnlimitedPopn`

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(100);

julia> u0_1 = simulationu0(; s=100_000, e=5, i_n=3, i_f=2);

julia> u0_2 = simulationu0(; s=200_000, e=5, i_n=3, i_f=2);

julia> gamma = 0.2; delta = 0.3; theta = 0.6; sigma = 0.5;

julia> _beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365);

julia> _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t);

julia> _beta2(t) = 0.72 * _beta1counter(t);

julia> s1 = packsimulationtuple( ;
       u0=u0_1, beta=_beta1, gamma, delta, theta, sigma, intervention=50,
       );

julia> s2 = packsimulationtuple( ;
       u0=u0_2, beta=_beta2, gamma, delta, theta, sigma, intervention=nothing,
       );

julia> packsimulationsunlimitedpopulation(rng, 100, s1, s2)
RenewalDiDDataUnlimitedPopn{Int64, InterventionMatrix{Int64}}
 observedcases:  [0 0; 0 0; … ; 6 1380; 9 1414]
 interventions:  [0 0; 0 0; … ; 1 0; 1 0] {duration 100, starttimes [50, nothing]}
 Ns:             unlimited
 exptdseedcases: [0.0 0.46548963316975533; 0.024913053640556182 0.0; … ; \
    0.4209971568176677 0.0; 0.5200181826119457 0.03451036683024468]
```
"""
function packsimulationsunlimitedpopulation(args...; kwargs...)
    return _packsimulations(_renewaldiddataunlimitedpopn_removens, args...; kwargs...)
end

function _renewaldiddataunlimitedpopn_removens(; Ns=nothing, kwargs...)
    return RenewalDiDDataUnlimitedPopn(; kwargs...)
end

function _packsimulations(f, duration, m1_args, args...; kwargs...)
    return _packsimulations(f, default_rng(), duration, m1_args, args...; kwargs...)
end

function _packsimulations(f, rng::AbstractRNG, duration, m1_args, args...; kwargs...)
    interventiontimes = Vector{Union{Int, Nothing}}(undef, 0)
    Ns = zeros(Int, 0)
    observedcases = zeros(Int, duration + 1, 0)
    interventiontimes, Ns, observedcases = _packsimulation!(
        rng, interventiontimes, Ns, observedcases, duration, m1_args, args...
    )
    interventions = InterventionMatrix{Int}(duration, interventiontimes)
    return f(; observedcases, interventions, Ns, kwargs...)
end

function _packsimulation!(rng, interventiontimes, Ns, observedcases, duration, m1_args)
    # simulation with one set of arguments
    return __packsimulation!(rng, interventiontimes, Ns, observedcases, duration, m1_args)
end

function _packsimulation!(
    rng, interventiontimes, Ns, observedcases, duration, m1_args, args...
)
    interventiontimes, Ns, observedcases = __packsimulation!(
        # use first set of arguments
        rng, interventiontimes, Ns, observedcases, duration, m1_args
    )
    # iterate until only one Tuple remains
    return _packsimulation!(rng, interventiontimes, Ns, observedcases, duration, args...)
end

function __packsimulation!(rng, interventiontimes, Ns, observedcases, duration, args::Tuple)
    u0, beta, gamma, delta, theta, sigma, intervention = args
    n = _n_seir(u0)
    cases = simulationcases(rng, duration, u0, beta, gamma, delta, theta, sigma)
    observedcases = hcat(observedcases, cases)
    push!(interventiontimes, intervention)
    push!(Ns, n)
    return (interventiontimes, Ns, observedcases)
end

"""
    packsimulationtuple(; u0, beta, gamma, delta, theta, sigma, intervention)

Return a tuple of arguments to be used with `packsimulations` or 
    `packsimulationsunlimitedpopulation`.

# Examples
```jldoctest
julia> u0 = simulationu0(; s=100_000, e=5, i_n=3, i_f=2);

julia> gamma = 0.2; delta = 0.3; theta = 0.6; sigma = 0.5;

julia> _beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365);

julia> _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t);

julia> packsimulationtuple( ;
       u0, beta=_beta1, gamma, delta, theta, sigma, intervention=50,
       )
([100000, 5, 3, 2, 0, 0, 0], _beta1, 0.2, 0.3, 0.6, 0.5, 50)
```
"""
function packsimulationtuple(; u0, beta, gamma, delta, theta, sigma, intervention)
    return (u0, beta, gamma, delta, theta, sigma, intervention)
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

_negativeuerror(u) = ArgumentError("$u, no members of u may be negative")
_negativedurationerror(d) = ArgumentError("$d: duration must be a positive integer")

function _parameterzeroerrorexception(x, t, symbol) 
    return ErrorException("$x: $symbol must be ≥ 0 (error at time $t)")
end

function _parametermaximumerrorexception(x, t, symbol, upper)
    return ErrorException("$x: $symbol must be ≤ $upper (error at time $t)")
end

function _recoveryerrorexception(timetodiagnosis, timetorecovery)
    m = "expected time to diagnosis ($timetodiagnosis) must be less than expected time to \
        recovery ($timetorecovery)"
    return ErrorException(m)
end

_simmatrixwidtherror(M) = ArgumentError("size $(size(M)): expecting a Matrix of width 7")

function _simulationu0_n_error(n, sm)
    m = "$n: when keyword argument `n` is supplied it must be at least as large as the sum \
        of the other compartment values provided ($sm)"
    return ArgumentError(m)
end

_simulationu0_negativeerror(x) = ArgumentError("$x: compartment sizes cannot be negative")

function _simulationu0_nr_error(n, calcn)
    m = "Inconsistent values of `n` and `r` supplied. Calculated n=$calcn but keyword \
        argument n=$n."
    return ArgumentError(m)
end

_ulengtherror(u) = ArgumentError("$u, u must be a vector of 7 integers")
