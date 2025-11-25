# generates simulated datasets

# `runsimulation` is the main function in this file. All parameters and compartments have
#  the meanings described in this function's documentation.
"""
    runsimulation([rng], duration, u0; beta, sigma, eta, phi)
    runsimulation([rng], duration, u0, beta, sigma, eta, phi)

Generate simulated data.

Model uses Gillespie's stochastic continuous time method with 6 compartments:
- `S`: susceptible 
- `E`: exposed (infected but not yet infectious)
- `I`: infectious and undiagnosed 
- `I′`: infectious and diagnosed 
- `R`: recovered

A 6th compartment in the model records the cumulative number of diagnoses.

Parameters `beta`, `sigma`, `eta` and `phi` may be numbers, or may be functions taking time 
as their only argument to give time-varying values. They can be entered as positional 
arguments or keyword arguments. All parameters must be non-negative.

# Arguments
- `rng::AbstractRNG=Random.default_rng()`: random number generator
- `duration::Integer`: duration of simulation
- `u0::AbstractVector{<:Integer}`: initial conditions, must be of length 7 describing
    compartments in the order above
- `beta`: transmission parameter
- `sigma`: rate of progression from exposed to infectious
- `eta`: recovery rate
- `phi`: proportion diagnosed; must satisfy `0 ≤ phi ≤ 1`

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
    # S   E   I   I′  R   cumulativediagnoses
     -1   1   0   0   0   0  # infections
      0  -1   1   0   0   0  # disease progression
      0   0  -1   1   0   1  # diagnosis
      0   0  -1   0   1   0  # recovery from I
      0   0   0  -1   1   0  # recovery from I′
]


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Simulation u0

"""
    simulationu0(; S, E, I, Iprime, R, N)

Produce vector of initial conditions for simulation.

# Keyword arguments 
The following represent the model compartments, 
- `S::Integer=0`: susceptible 
- `E::Integer=0`: exposed (infected but not yet infectious)
- `I::Integer=0`: infectious and undiagnosed 
- `Iprime::Integer=0`: infectious and diagnosed 
- `R::Union{<:Integer, Automatic}=automatic`: recovered

The argument `N::Union{<:Integer, Nothing}=nothing` is the population size. If `N` is 
    provided and `R` is not then the remaining population not assigned to other compartments 
    will be assigned to `R`. An error is thrown if a unique non-negative value of `R` cannot
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

julia> simulationu0(; S=99, E=1)
7-element Vector{Int64}:
 99
  1
  0
  0
  0
  0

julia> simulationu0(; S=99, E=1, N=200)
7-element Vector{Int64}:
  99
   1
   0
   0
 100
   0

julia> simulationu0(; S=99, E=1, R=100, N=250)
ERROR: ArgumentError: Inconsistent values of `n` and `r` supplied. Calculated n=200 but \
    keyword argument n=250.
```
"""
function simulationu0( ; 
    S::Integer=0, 
    E::Integer=0, 
    I::Integer=0, 
    Iprime::Integer=0, 
    R::Union{<:Integer, Automatic}=automatic,
    N::Union{<:Integer, Nothing}=nothing,
)
    return _simulationu0(S, E, I, Iprime, R, N)
end

function _simulationu0(S, E, I, Iprime, ::Automatic, ::Nothing)
    return _simulationu0(S, E, I, Iprime, 0, nothing)
end

function _simulationu0(S, E, I, Iprime, ::Automatic, N::Integer)
    R = N - (S + E + I + Iprime)
    R >= 0 || throw(_simulationu0_n_error(N, N - R))
    return _simulationu0(S, E, I, Iprime, R, nothing)
end

function _simulationu0(S, E, I, Iprime, R::Integer, N::Integer)
    S + E + I + Iprime + R == N || throw(_simulationu0_nr_error(N, S + E + I + Iprime + R))
    return _simulationu0(S, E, I, Iprime, R, nothing)
end

function _simulationu0(S, E, I, Iprime, R::Integer, ::Nothing)
    u0 = [S, E, I, Iprime, R, 0]
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
    x < upper || throw(_parametermaximumerrorexception(x, t, symbol, upper))
    return x
end

_pbeta(beta, t) = _parameter(beta, t, :beta)
_peta(eta, t) = _parameter(eta, t, :eta)
_psigma(sigma, t) = _parameter(sigma, t, :sigma)
_pphi(phi, t) = _parameter(phi, t, :phi; upper=1)  # phi is a proportion 0 < φ < 1

## Force of infection 

# i_n, i_f and i_d assumed to be equally infectious
_foi(beta, I, Iprime, N, t) = _pbeta(beta, t) * (I + Iprime) / N

## Event rates

_simulatedinfections(beta, S, I, Iprime, N, t) = S * _foi(beta, I, Iprime, N, t)
_diseaseprogression(E, sigma, t) = E * _psigma(sigma, t)
_diagnosis(I, eta, phi, t) = I * _peta(eta, t) * _pphi(phi, t) / (1 - _pphi(phi, t))
_recovery(x, eta, t) = x * _peta(eta, t)

function _n_seir(u::AbstractVector{<:Integer})
    length(u) == 6 || throw(_ulengtherror(u))
    minimum(u) >= 0 || throw(_negativeuerror(u))
    return sum(@view u[1:5])  # S, E, I, Iprime, R
end

function _seirrates(u::AbstractVector{<:Integer}, t, beta, sigma, eta, phi)
    N = _n_seir(u)
    S, E, I, Iprime, R, = u
    return [
        _simulatedinfections(beta, S, I, Iprime, N, t),  # infections
        _diseaseprogression(E, sigma, t),  # disease progression 
        _diagnosis(I, eta, phi, t),  # diagnosis
        _recovery(I, eta, t),  # recovery from I
        _recovery(Iprime, eta, t),  # recovery from I′
    ]
end

## Effect the next event

_tstep(rng, rates) = -log(rand(rng)) / sum(rates)
_nextevent(rates) = _nextevent(default_rng(), rates)  
# `_nextevent(rates)` is not used by any functions except the tests -- may be worth changing 
# the tests and removing this method
_nextevent(rng, rates) = sample(rng, eachindex(rates), Weights(rates))
_updateevent!(u, nextevent) = u .+= _SEIREVENTSMATRIX[nextevent, :]

## Simulate a day

_simulateday!(args...) = _simulateday!(default_rng(), args...)

function _simulateday!(rng::AbstractRNG, u, t, beta, sigma, eta, phi)
    nextday = t + 1 

    while t < nextday
        rates = _seirrates(u, t, beta, sigma, eta, phi)
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

function runsimulation(rng::AbstractRNG, duration, u0; beta, sigma, eta, phi)
    return _runsimulation(rng, duration, u0, beta, sigma, eta, phi)
end

function runsimulation(rng::AbstractRNG, duration, u0, beta, sigma, eta, phi)
    return _runsimulation(rng, duration, u0, beta, sigma, eta, phi) 
end

function _runsimulation(rng, duration, u0, beta, sigma, eta, phi)
    duration >= 1 || throw(_negativedurationerror(duration))
    output = zeros(Int, duration + 1, 6)
    u = deepcopy(u0)
    output[1, :] = u
    _runsimulationdays!(rng, output, duration, u, beta, sigma, eta, phi)
    return output 
end

function _runsimulationdays!(rng, output, duration, u, beta, sigma, eta, phi)
    for t in 1:duration 
        _simulateday!(rng, u, t, beta, sigma, eta, phi)
        output[t+1, :] = u
    end
    return nothing
end

"""
    simulationcases(M::Matrix)
    simulationcases([rng], duration, u0, beta, sigma, eta, phi)
    simulationcases([rng], duration, u0; beta, sigma, eta, phi)

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
    size(M, 2) == 6 || throw(_simmatrixwidtherror(M))
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
    xs[1] == M[1, 6]
    for t in 2:(duration + 1)
        xs[t] = M[t, 6] - M[t-1, 6]
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
    return __packsimulationsinterventiontimesarray(m1_args[6])
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
    u0, beta, sigma, eta, phi, intervention = args
    N = _n_seir(u0)
    cases = simulationcases(rng, duration, u0, beta, sigma, eta, phi)
    observedcases = hcat(observedcases, cases)
    interventiontimes = _packsimulationinterventiontimes(interventiontimes, intervention)
    push!(Ns, N)
    return (interventiontimes, Ns, observedcases)
end

function _packsimulationinterventiontimes(interventiontimes::Matrix, intervention::Vector)
    return vcat(interventiontimes, permutedims(intervention))
end

function _packsimulationinterventiontimes(interventiontimes::Vector, intervention)
    return vcat(interventiontimes, intervention)
end

const _SIMARGSORDER = [:u0, :beta, :eta, :sigma, :phi]

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
function packsimulationtuple(; u0, beta, sigma, eta, phi, intervention)
    return (u0, beta, sigma, eta, phi, intervention)
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

_negativeuerror(u) = ArgumentError("$u, no members of u may be negative")
_negativedurationerror(d) = ArgumentError("$d: duration must be a positive integer")

function _parameterzeroerrorexception(x, t, symbol) 
    return ErrorException("$x: $symbol must be ≥ 0 (error at time $t)")
end

function _parametermaximumerrorexception(x, t, symbol, upper)
    return ErrorException("$x: $symbol must be < $upper (error at time $t)")
end

function _recoveryerrorexception(timetodiagnosis, timetorecovery)
    m = "expected time to diagnosis ($timetodiagnosis) must be less than expected time to \
        recovery ($timetorecovery)"
    return ErrorException(m)
end

_simmatrixwidtherror(M) = ArgumentError("size $(size(M)): expecting a Matrix of width 6")

function _simulationu0_n_error(n, sm)
    m = "$n: when keyword argument `n` is supplied it must be at least as large as the sum \
        of the other compartment values provided ($sm)"
    return ArgumentError(m)
end

_simulationu0_negativeerror(x) = ArgumentError("$x: compartment sizes cannot be negative")

function _simulationu0_nr_error(N, calcn)
    msg = "Inconsistent values of `N` and `R` supplied. Calculated N=$calcn but keyword \
        argument N=$N."
    return ArgumentError(msg)
end

_ulengtherror(u) = ArgumentError("$u, u must be a vector of 6 integers")
