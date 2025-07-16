const _SEIREVENTSMATRIX = [
    # s   e   i   i′  r   cumulativediagnoses
     -1   1   0   0   0   0  # infections
      0  -1   1   0   0   0  # disease progression
      0   0  -1   1   0   1  # diagnosis
      0   0  -1   0   1   0  # recovery from i
      0   0   0  -1   1   0  # recovery from i′
]

function _parameter(x::Number, ::Any, symbol)
    x >= 0 || throw(ErrorException("$x: $symbol must be ≥ 0"))
    return x
end

_parameter(f::Function, t, symbol) = _parameter(f(t), nothing, symbol)
_parameter(::Nothing, ::Any, ::Any) = 0
_parameter(x, t) = _parameter(x, t, :parameter)
_parametermultiplication(x, args...) = x * _parameter(args...)
_foi(beta, i, i_prime, n, t) = _parameter(beta, t, :β) * (i + i_prime) / n
_simulatedinfections(beta, s, i, i_prime, n, t) = s * _foi(beta, i, i_prime, n, t)
_diseaseprogression(x, sigma, t) = _parametermultiplication(x, sigma, t, :ς)
_recovery(x, gamma, t) = _parametermultiplication(x, gamma, t, :ς)

function _diagnosis(x, gamma, theta, t)
    θ = _parameter(theta, t, :θ) 
    θ < 1 || throw(ArgumentError("$θ, theta must be strictly less than 1"))
    return θ * _parameter(gamma, t, :γ) * x / (1 - θ)
end

function _seirrates(u::AbstractVector{<:Integer}, t, beta, sigma, gamma, theta)
    length(u) == 6 || throw(ArgumentError("$u, u must be a vector of 6 integers"))
    n = sum(@view u[1:5])  # s, e, i, i′, r
    s, e, i, i_prime, = u
    return [
        _simulatedinfections(beta, s, i, i_prime, n, t),  # leave s enter e
        _diseaseprogression(e, sigma, t),  # leave e and enter i
        _diagnosis(i, gamma, theta, t),  # leave i and enter i′
        _recovery(i, gamma, t),  # leave i and enter r
        _recovery(i_prime, gamma, t)  # leave i_prime and enter r
    ]
end

_tstep(rng, rates) = -log(rand(rng)) / sum(rates)

_nexteventtime(t, rates) = _nexteventtime(default_rng(), t, rates)

function _nexteventtime(rng::AbstractRNG, t, rates)
    tstep = _tstep(rng, rates)
    return t + tstep
end

_nextevent(rates) = _nextevent(default_rng(), rates)
_nextevent(rng::AbstractRNG, rates) = sample(rng, eachindex(rates), Weights(rates))
_updateevent!(u, nextevent) = u .+= _SEIREVENTSMATRIX[nextevent, :]

function _simulateday!(u::AbstractVector{<:Integer}, t, beta, sigma, gamma, theta)
    return _simulateday!(default_rng(), u, t, beta, sigma, gamma, theta)
end

function _simulateday!(
    rng::AbstractRNG, u::AbstractVector{<:Integer}, t, beta, sigma, gamma, theta
)
    nextday = t + 1 
    while t < nextday
        rates = _seirrates(u, t, beta, sigma, gamma, theta)
        nextevent = _nextevent(rng, rates)
        t += _tstep(rng, rates)
        if t < nextday 
            _updateevent!(u, nextevent) 
        end
    end
    return u 
end

"""
    runsimulation([rng], duration, u0, beta, sigma, gamma, theta)

Runs simulation for 1 day.

Parameters `beta`, `sigma`, `gamma` and `theta` may be numbers or may be time-varying 
functions taking a single argument `t`. `u` is mutated to give the conditions at the end of 
the day. 

# Arguments
- `rng::AbstractRNG=Random.default_rng()`: random number generator
- `duration::Integer`: duration of simulation 
- `u::AbstractVector{<:Integer}`: current conditions, must be of length 6 (susceptible, 
    exposed, infectious [undiagnosed], infectious [diagnosed], recovered, cumulative 
    diagnoses)
- `beta`: transmission parameter
- `sigma`: rate of progression from exposed to infectious
- `gamma`: recovery rate
- `theta`: proportion diagnosed (must satisfy `0 < θ ≤ 1`)

# Returns
Returns a matrix of length `duration + 1` giving numbers in each compartment at the end of 
each day. The first row is the conditions in `u0`.
"""
function runsimulation(duration, u0::AbstractVector{<:Integer}, beta, sigma, gamma, theta)
    return runsimulation(default_rng(), duration, u0, beta, sigma, gamma, theta) 
end

function runsimulation(
    rng::AbstractRNG, duration, u0::AbstractVector{<:Integer}, beta, sigma, gamma, theta
)
    length(u0) == 6 || throw(ArgumentError("$u0, u0 must be a vector of 6 integers"))
    output = zeros(Int, duration + 1, 6)
    u = deepcopy(u0)
    output[1, :] = u
    for t in 1:duration 
        _simulateday!(rng, u, t, beta, sigma, gamma, theta)
        output[t+1, :] = u
    end
    return output 
end

function simulationcases(M::Matrix)
    size(M, 2) == 6 || throw(ArgumentError("size $(size(M)): expecting a Matrix of width 6"))
    duration = size(M, 1) - 1
    return _simulationcases(M, duration)
end

function simulationcases(duration, args...)
    M = runsimulation(duration, args...)
    return _simulationcases(M, duration)
end

function simulationcases(rng::AbstractRNG, duration, args...)
    M = runsimulation(rng, duration, args...)
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

# not yet tested 
function packsimulations(duration, m1_args, args...)
    return packsimulations(default_rng(), duration, m1_args, args...)
end

function packsimulations(rng::AbstractRNG, duration, m1_args, args...)
    interventiontimes = Vector{Union{Int, Nothing}}(undef, 0)
    Ns = zeros(Int, 0)
    observedcases = zeros(Int, duration + 1, 0)
    interventiontimes, Ns, observedcases = _packsimulation!(
        interventiontimes, Ns, rng, observedcases, duration, m1_args, args...
    )

    return Dict(
        :observedcases => observedcases,
        :interventions => InterventionMatrix{Int}(duration, interventiontimes),
        :Ns => Ns,
    )
end

function _packsimulation!(interventiontimes, Ns, rng, observedcases, duration, m1_args)
    return __packsimulation!(interventiontimes, Ns, rng, observedcases, duration, m1_args)
end

function _packsimulation!(
    interventiontimes, Ns, rng, observedcases, duration, m1_args, args...
)
    interventiontimes, Ns, observedcases = __packsimulation!(
        interventiontimes, Ns, rng, observedcases, duration, m1_args
    )
    return _packsimulation!(interventiontimes, Ns, rng, observedcases, duration, args...)
end

function __packsimulation!(
    interventiontimes, Ns, rng::AbstractRNG, observedcases, duration, args::Tuple
)
    u0, beta, sigma, gamma, theta, intervention = args
    n = sum(@view u0[1:5])
    cases = simulationcases(rng, duration, u0, beta, sigma, gamma, theta)
    observedcases = hcat(observedcases, cases)
    push!(interventiontimes, intervention)
    push!(Ns, n)
    return (interventiontimes, Ns, observedcases)
end
