# generates simulated datasets

# `runsimulation` is the main function in this file. All parameters and compartments have
#  the meanings described in this function's documentation.
"""
    runsimulation([rng], duration, u0; beta, mu, delta, psi, kappa)
    runsimulation([rng], duration, u0, beta, mu, delta, psi, kappa)

Generate simulated data.

Model uses Gillespie's stochastic continuous time method with 6 compartments:
- `s`: susceptible 
- `e`: exposed (infected but not yet infectious)
- `i_n`: infectious and will never be diagnosed 
- `i_f`: infectious and will be diagnosed before recovery 
- `i_d`: infectious and diagnosed 
- `r`: recovered

A 7th compartment in the model records the cumulative number of diagnoses.

Parameters `beta`, `mu`, `delta` `psi` and `kappa` may be numbers, or may be functions 
taking time as their only argument to give time-varying values. They can be entered as 
positional arguments or keyword arguments. All parameters must be non-negative.

# Arguments
- `rng::AbstractRNG=Random.default_rng()`: random number generator
- `duration::Integer`: duration of simulation
- `u0::AbstractVector{<:Integer}`: initial conditions, must be of length 7 describing
    compartments in the order above
- `beta`: transmission parameter
- `mu`: recovery rate
- `delta`: diagnosis rate (i.e. for those who will be diagnosed, how quickly does it 
    happen?); must satisfy `delta > mu` (i.e. must be diagnosed before the end of the 
    infectious period)
- `psi`: proportion diagnosed; must satisfy `0 ≤ psi ≤ 1`
- `kappa`: rate of progression from exposed to infectious

# Returns
Returns a matrix of height `duration + 1` giving numbers in each compartment at the end of 
each day. The first row is the conditions supplied in `u0`.

See `packsimulations` to generate a `SimulationData` struct.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(1);

julia> u0 = simulationu0(; s=100, e=5, n=200);

julia> runsimulation(rng, 10, u0; beta=0.6, mu=0.25, delta=0.33, psi=0.5, kappa=0.4)
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

julia> runsimulation(rng, 10, u0; beta=mybetafunc, mu=0.25, delta=0.33, psi=0.5, kappa=0.4)
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
julia> simulationu0()
7-element Vector{Int64}:
 0
 0
 0
 0
 0
 0
 0

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
_pmu(mu, t) = _parameter(mu, t, :mu)
_pkappa(kappa, t) = _parameter(kappa, t, :kappa)
_ppsi(psi, t) = _parameter(psi, t, :psi; upper=1)  # psi is a proportion 0 ≤ ψ ≤ 1

## Force of infection 

# i_n, i_f and i_d assumed to be equally infectious
_foi(beta, i_n, i_f, i_d, n, t) = _pbeta(beta, t) * (i_n + i_f + i_d) / n

## Event rates

_simulatedinfections(beta, s, i_n, i_f, i_d, n, t) = s * _foi(beta, i_n, i_f, i_d, n, t)

function _diseaseprogression_to_i_n(x, psi, kappa, t)
    return x * (1 - _ppsi(psi, t)) * _pkappa(kappa, t)
end

_diseaseprogression_to_i_f(x, psi, kappa, t) = x * _ppsi(psi, t) * _pkappa(kappa, t)
_diagnosis(x, delta, t) = x * _pdelta(delta, t)
_recovery(x, mu, t) = x * _pmu(mu, t)

# diagnosed individuals have already been infectious for a period represented by `1/delta`
# so recover after a period `1/mu - 1/delta`
function _recovery(x, mu, delta, t) 
    timetodiagnosis = 1 / _pdelta(delta, t)
    timetorecovery = 1 / _pmu(mu, t)
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

function _seirrates(u::AbstractVector{<:Integer}, t, beta, mu, delta, psi, kappa)
    n = _n_seir(u)
    s, e, i_n, i_f, i_d, = u
    return [
        _simulatedinfections(beta, s, i_n, i_f, i_d, n, t),  # infections
        _diseaseprogression_to_i_n(e, psi, kappa, t),  # disease progression, not diagnosed
        _diseaseprogression_to_i_f(e, psi, kappa, t),  # disease progression, to be diagnosed
        _diagnosis(i_f, delta, t),  # diagnosis
        _recovery(i_n, mu, t),  # recovery from i_n
        _recovery(i_d, mu, delta, t)  # recovery from i_d
    ]
end

## Effect the next event

_tstep(rng, rates) = -log(rand(rng)) / sum(rates)
_nextevent(rates) = _nextevent(default_rng(), rates)  
# `_nextevent(rates)` is not used by any functions expect the tests -- may be worth changing 
# the tests and removing this method
_nextevent(rng, rates) = sample(rng, eachindex(rates), Weights(rates))
_updateevent!(u, nextevent) = u .+= _SEIREVENTSMATRIX[nextevent, :]

## Simulate a day

_simulateday!(args...) = _simulateday!(default_rng(), args...)

function _simulateday!(rng::AbstractRNG, u, t, beta, mu, delta, psi, kappa)
    nextday = t + 1 
    while t < nextday
        rates = _seirrates(u, t, beta, mu, delta, psi, kappa)
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

function runsimulation(rng::AbstractRNG, duration, u0; beta, mu, delta, psi, kappa)
    return _runsimulation(rng, duration, u0, beta, mu, delta, psi, kappa)
end

function runsimulation(rng::AbstractRNG, duration, u0, beta, mu, delta, psi, kappa)
    return _runsimulation(rng, duration, u0, beta, mu, delta, psi, kappa) 
end

function _runsimulation(rng, duration, u0, beta, mu, delta, psi, kappa)
    duration >= 1 || throw(_negativedurationerror(duration))
    output = zeros(Int, duration + 1, 7)
    u = deepcopy(u0)
    output[1, :] = u
    _runsimulationdays!(rng, output, duration, u, beta, mu, delta, psi, kappa)
    return output 
end

function _runsimulationdays!(rng, output, duration, u, beta, mu, delta, psi, kappa)
    for t in 1:duration 
        _simulateday!(rng, u, t, beta, mu, delta, psi, kappa)
        output[t+1, :] = u
    end
    return nothing
end

"""
    simulationcases(M::Matrix)
    simulationcases([rng], duration, u0, beta, mu, delta, psi, kappa)
    simulationcases([rng], duration, u0; beta, mu, delta, psi, kappa)

Produce a vector of simulated diagnosed infections.

Either takes a matrix that is the output of `runsimulation` or a set of arguments that gets 
passed to `runsimulation`. 

See `runsimulation` for more details.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng1 = StableRNG(1);

julia> u0 = simulationu0(; s=100, e=5, n=200);

julia> c = simulationcases(rng1, 10, u0; beta=0.6, mu=0.25, delta=0.33, psi=0.5, kappa=0.4)
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

julia> rng2 = StableRNG(1);

julia> sim = runsimulation(rng2, 10, u0; beta=0.6, mu=0.25, delta=0.33, psi=0.5, kappa=0.4);

julia> simulationcases(sim) == c
true
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

Run simulations and collate results into a `SimulationData` struct.

# Arguments
- `rng::AbstractRNG=Random.default_rng()``: random number generator
- `duration`: duration of the simulation.
- `m1_args`: a `Tuple` containing `{u0, beta, mu, delta, psi, kappa, intervention}` in 
    order for the first group's simulation (the function `packsimulationtuple` can be used 
    to generate this with keyword arguments)
- `args...`: equivalent Tuples for remaining groups
- keyword arguments are passed to `SimulationData`

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(100);

julia> u0_1 = simulationu0(; s=100_000, e=5, i_n=3, i_f=2);

julia> u0_2 = simulationu0(; s=200_000, e=5, i_n=3, i_f=2);

julia> mu = 0.2; delta = 0.3; psi = 0.6; kappa = 0.5;

julia> beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365);  # intervention-free counterfactual

julia> beta1(t) = t >= 50 ? 0.8 * beta1counter(t) : beta1counter(t);

julia> beta2(t) = 0.72 * beta1counter(t);

julia> s1 = packsimulationtuple(; u0=u0_1, beta=beta1, mu, delta, psi, kappa, intervention=50);

julia> s2 = packsimulationtuple(; u0=u0_2, beta=beta2, mu, delta, psi, kappa, intervention=nothing);

julia> packsimulations(rng, 100, s1, s2)
SimulationData{Int64, InterventionMatrix{Int64}, Vector{Int64}}
 observedcases:  [0 0; 0 0; … ; 6 1380; 9 1414]
 interventions:  [0 0; 0 0; … ; 1 0; 1 0] {duration 100, starttimes [50, nothing]}
 Ns:             [100010, 200010]
 exptdseedcases: [0.0 0.0; 0.024913053640556182 0.0; … ; 0.4209971568176677 0.0; 0.5200181826119457 0.5]
```
"""
function packsimulations(args...; kwargs...)
    return _packsimulations(args...; kwargs...)
end

function _packsimulations(duration, m1_args, args...; kwargs...)
    return _packsimulations(default_rng(), duration, m1_args, args...; kwargs...)
end

function _packsimulations(
    rng::AbstractRNG, duration, m1_args, args...; 
    unlimitedpop=false, kwargs...
)
    interventiontimes = _packsimulationsinterventiontimesarray(m1_args)
    Ns = zeros(Int, 0)
    observedcases = zeros(Int, duration + 1, 0)
    interventiontimes, Ns, observedcases = _packsimulation!(
        rng, interventiontimes, Ns, observedcases, duration, m1_args, args...
    )
    interventions = _siminterventionarray(duration, interventiontimes)
    u0s = _simargs(:u0, m1_args, args...)
    if unlimitedpop 
        return SimulationData(; observedcases, interventions, Ns=nothing, rng, u0s, kwargs...)
    else 
        return SimulationData(; observedcases, interventions, Ns, rng, u0s, kwargs...)
    end
end

function _siminterventionarray(duration, interventiontimes::Vector)
    return InterventionMatrix{Int}(duration, interventiontimes)
end

function _siminterventionarray(duration, interventiontimes::Matrix)
    return InterventionArray{Int}(duration, interventiontimes)
end

function _packsimulationsinterventiontimesarray(m1_args)
    return __packsimulationsinterventiontimesarray(m1_args[7])
end

function __packsimulationsinterventiontimesarray(::Number)
    return __packsimulationsinterventiontimesarray(nothing)
end

__packsimulationsinterventiontimesarray(::Nothing) = Vector{Union{Int, Nothing}}(undef, 0)

function __packsimulationsinterventiontimesarray(v::Vector)
    return Matrix{Union{Int, Nothing}}(undef, 0, length(v))
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
    u0, beta, mu, delta, psi, kappa, intervention = args
    n = _n_seir(u0)
    cases = simulationcases(rng, duration, u0, beta, mu, delta, psi, kappa)
    observedcases = hcat(observedcases, cases)
    interventiontimes = _packsimulationinterventiontimes(interventiontimes, intervention)
    push!(Ns, n)
    return (interventiontimes, Ns, observedcases)
end

function _packsimulationinterventiontimes(interventiontimes::Matrix, intervention::Vector)
    return vcat(interventiontimes, permutedims(intervention))
end

function _packsimulationinterventiontimes(interventiontimes::Vector, intervention)
    return vcat(interventiontimes, intervention)
end

const _SIMARGSORDER = [:u0, :beta, :mu, :delta, :psi, :kappa]

function _simargs(parameter::Symbol, args...)
    index = findfirst(x -> x == parameter, _SIMARGSORDER)
    return _simargs(index::Int, args...)
end

_simargs(index::Int, args...) = [a[index] for a in args]

"""
    packsimulationtuple(; u0, beta, mu, delta, psi, kappa, intervention)

Return a tuple of arguments to be used with `packsimulations`.

# Examples
```jldoctest
julia> u0 = simulationu0(; s=100_000, e=5, i_n=3, i_f=2);

julia> beta = 0.4; mu = 0.2; delta = 0.3; psi = 0.6; kappa = 0.5;

julia> packsimulationtuple(; u0, beta, mu, delta, psi, kappa, intervention=50)
([100000, 5, 3, 2, 0, 0, 0], 0.4, 0.2, 0.3, 0.6, 0.5, 50)
```
"""
function packsimulationtuple(; u0, beta, mu, delta, psi, kappa, intervention)
    return (u0, beta, mu, delta, psi, kappa, intervention)
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
