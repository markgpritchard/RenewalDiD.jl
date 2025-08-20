# data types used by `renewaldid`

"""
    AbstractRenewalDiDData{S, T}

Abstract type of containers for data passed to `renewaldid`.

Two subtypes are defined:
- `RenewalDiDData{S, T}`, which allows tracking of numbers susceptible
- `RenewalDiDDataUnlimitedPopn{S, T}`, which assumes an unlimited population that is not
    depleted by infections
"""
abstract type AbstractRenewalDiDData{S, T} end

"""
    RenewalDiDData{S, T}

Container for data passed to `renewaldid`.

# Fields
- `observedcases::Matrix{S}`: matrix of observed cases, each group is a separate column
- `interventions::T`: array of intervention times. Time changes by row, group by column, and
    intervention in dimension 3 (see `InterventionMatrix` and `InterventionArray`) 
- `Ns::Vector{Int}`: population size of each group
- `exptdseedcases::Matrix{Float64}`: matrix of infections up to time `t=0` that seeds
    subsequent infection events

# Constructors

    RenewalDiDData(; observedcases, interventions, Ns, exptdseedcases=nothing, <keyword \
        arguments>)

The constructor takes all arguments as keyword arguments.  

The matrix `exptdseedcases` may be supplied or one is generated automatically. Remaining 
    keyword arguments are passed to `expectedseedcases`, with that function taking its 
    default arguments if not specified.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(10);

julia> observedcases = rand(rng, Poisson(10), 11, 3);  # 3 groups for times 0 to 10

julia> interventions = InterventionMatrix(10, 3, 7, nothing);

julia> Ns = 1000 .* ones(Int, 3);

julia> RenewalDiDData(; observedcases, interventions, Ns)
RenewalDiDData{Int64, InterventionMatrix{Int64}}
 observedcases:  [10 8 15; 20 8 11; … ; 10 6 14; 9 9 8]
 interventions:  [0 0 0; 0 0 0; … ; 1 1 0; 1 1 0] {duration 10, starttimes [3, 7, nothing]}
 Ns:             [1000, 1000, 1000]
 exptdseedcases: [1.6236225474260566 1.4880770554298328 1.8607523407150066; \
    1.7226435732203347 1.587098081224111 1.9597733665092847; … ; 2.1187276763974463 \
    1.9831821844012225 2.3558574696863963; 2.217748702191724 2.0822032101955004 \
    2.454878495480674]
```
""" 
@auto_hash_equals struct RenewalDiDData{S, T} <: AbstractRenewalDiDData{S, T}
    observedcases::Matrix{S}
    interventions::T
    Ns::Vector{Int}
    exptdseedcases::Matrix{Float64}
    id::String

    function RenewalDiDData(
        observedcases::Matrix{S}, 
        interventions::T, 
        Ns::Vector{<:Number}, 
        exptdseedcases::Matrix{<:Number},
        id::String=""
    ) where {S <:Number, T <:AbstractArray}
        _diddataassertions(observedcases, interventions, Ns, exptdseedcases)
        return new{S, T}(
            observedcases, 
            interventions, 
            convert.(Int, Ns), 
            convert.(Float64, exptdseedcases),
            id
        )
    end
end

function RenewalDiDData( ; 
    observedcases, 
    interventions, 
    Ns, 
    exptdseedcases=nothing,
    n_seeds=DEFAULT_SEEDMATRIX_HEIGHT, 
    doubletime=automatic, 
    sampletime=automatic, 
    minvalue=DEFAULT_SEEDMATRIX_MINVALUE,
    id="",
)
    newexptdseedcases = _expectedseedcasesifneeded(
        exptdseedcases, observedcases, n_seeds; 
        doubletime, sampletime, minvalue
    )
    return RenewalDiDData(observedcases, interventions, Ns, newexptdseedcases, id)
end

"""
    RenewalDiDDataUnlimitedPopn{S, T}

Container for data passed to `renewaldid`.

Differs from `RenewalDiDData` in that it tells `renewaldid` to assume an inexhaustible 
    population so changes in susceptibility are not explicitly tracked.

# Fields
- `observedcases::Matrix{S}`: matrix of observed cases, each group is a separate column
- `interventions::T`: matrix of intervention times, each group is a separate column 
- `exptdseedcases::Matrix{Float64}`: matrix of infections up to time `t=0` that seeds
     subsequent infection events

# Constructors

    RenewalDiDDataUnlimitedPopn(; observedcases, interventions, exptdseedcases=nothing, \
        <keyword arguments>)

The constructor takes all arguments as keyword arguments.  

The matrix `exptdseedcases` may be supplied or one is generated automatically. Remaining 
    keyword arguments are passed to `expectedseedcases`, with that function taking its 
    default arguments if not specified.

If a keyword argument `Ns` is supplied it is not used and a warning is generated.

# Examples
```jldoctest
julia> using StableRNGs

julia> rng = StableRNG(10);

julia> observedcases = rand(rng, Poisson(10), 11, 3);  # 3 groups for times 0 to 10

julia> interventions = InterventionMatrix(10, 3, 7, nothing);

julia> RenewalDiDDataUnlimitedPopn(; observedcases, interventions)
RenewalDiDDataUnlimitedPopn{Int64, InterventionMatrix{Int64}}
 observedcases:  [10 8 15; 20 8 11; … ; 10 6 14; 9 9 8]
 interventions:  [0 0 0; 0 0 0; … ; 1 1 0; 1 1 0] {duration 10, starttimes [3, 7, nothing]}
 Ns:             unlimited
 exptdseedcases: [1.6236225474260566 1.4880770554298328 1.8607523407150066; \
    1.7226435732203347 1.587098081224111 1.9597733665092847; … ; 2.1187276763974463 \
    1.9831821844012225 2.3558574696863963; 2.217748702191724 2.0822032101955004 2.454878495480674]

julia> v = RenewalDiDDataUnlimitedPopn(; observedcases, interventions, Ns=(1000 .* ones(Int, 3)));
┌ Warning: [1000, 1000, 1000]: keyword argument `Ns` not used in constructing `RenewalDiDDataUnlimitedPopn`
└ @ RenewalDiD 

julia> isnothing(v.Ns)
true
```
""" 
@auto_hash_equals struct RenewalDiDDataUnlimitedPopn{S, T} <: AbstractRenewalDiDData{S, T}
    observedcases::Matrix{S}
    interventions::T
    exptdseedcases::Matrix{Float64}
    id::String

    function RenewalDiDDataUnlimitedPopn(
        observedcases::Matrix{S}, 
        interventions::T, 
        exptdseedcases::Matrix{<:Number},
        id::String=""
    ) where {S <:Number, T <:AbstractArray}
        _diddataassertions(observedcases, interventions, exptdseedcases)
        return new{S, T}(
            observedcases, 
            interventions, 
            convert.(Float64, exptdseedcases),
            id
        )
    end
end

function RenewalDiDDataUnlimitedPopn( ; 
    observedcases, 
    interventions, 
    Ns=nothing,
    exptdseedcases=nothing,
    n_seeds=DEFAULT_SEEDMATRIX_HEIGHT, 
    doubletime=automatic, 
    sampletime=automatic, 
    minvalue=DEFAULT_SEEDMATRIX_MINVALUE,
    id="",
)
    return _RenewalDiDDataUnlimitedPopn(
        observedcases, 
        interventions, 
        Ns,
        exptdseedcases,
        id,
        n_seeds, 
        doubletime, 
        sampletime, 
        minvalue,
    )
end

function _RenewalDiDDataUnlimitedPopn(
    observedcases, 
    interventions, 
    Ns,
    exptdseedcases,
    id,
    n_seeds, 
    doubletime, 
    sampletime, 
    minvalue,
)
    @warn "$Ns: keyword argument `Ns` not used in constructing `RenewalDiDDataUnlimitedPopn`"
    return _RenewalDiDDataUnlimitedPopn(
        observedcases, 
        interventions, 
        nothing,
        exptdseedcases,
        id,
        n_seeds, 
        doubletime, 
        sampletime, 
        minvalue,
    )
end

function _RenewalDiDDataUnlimitedPopn(
    observedcases, 
    interventions, 
    ::Nothing,
    exptdseedcases,
    id,
    n_seeds, 
    doubletime, 
    sampletime, 
    minvalue,
)
    newexptdseedcases = _expectedseedcasesifneeded(
        exptdseedcases, observedcases, n_seeds; 
        doubletime, sampletime, minvalue
    )
    return RenewalDiDDataUnlimitedPopn(observedcases, interventions, newexptdseedcases, id)
end

## simulation data

@auto_hash_equals struct SimulationData{S, T, U, V, W, X, Y, Z} <: AbstractRenewalDiDData{S, T}
    observedcases::Matrix{S}
    interventions::T
    Ns::Vector{Int}
    exptdseedcases::Matrix{Float64}
    rng::U
    u0s::Vector{Vector{Int}}
    betas::Vector{V}
    mus::Vector{W} 
    deltas::Vector{X}
    psis::Vector{Y}
    kappas::Vector{Z}
    id::String

    function SimulationData(
        observedcases::Matrix{S}, 
        interventions::T, 
        Ns::Vector{Int}, 
        exptdseedcases::Matrix{<:Number}, 
        rng::U, 
        u0s::Vector{Vector{Int}}, 
        betas::Vector{V}, 
        mus::Vector{W} , 
        deltas::Vector{X}, 
        psis::Vector{Y}, 
        kappas::Vector{Z}, 
        id::String=""
    ) where {
        S <: Number, 
        T <: AbstractArray, 
        U <: AbstractRNG, 
        V <: Any, 
        W <: Any, 
        X <: Any, 
        Y <: Any, 
        Z <: Any
    }
        _diddataassertions(observedcases, interventions, Ns, exptdseedcases)
        return new{S, T, U, V, W, X, Y, Z}(
            observedcases, 
            interventions, 
            Ns, 
            convert.(Float64, exptdseedcases),
            rng,
            u0s,
            betas,
            mus,
            deltas,
            psis,
            kappas,
            id
        )
    end
end

function SimulationData( ; 
    observedcases, 
    interventions, 
    Ns, 
    rng, 
    u0s, 
    betas, 
    mus, 
    deltas, 
    psis, 
    kappas,
    exptdseedcases=nothing,
    n_seeds=DEFAULT_SEEDMATRIX_HEIGHT, 
    doubletime=automatic, 
    sampletime=automatic, 
    minvalue=DEFAULT_SEEDMATRIX_MINVALUE,
    id="",
)
    newexptdseedcases = _expectedseedcasesifneeded(
        exptdseedcases, observedcases, n_seeds; 
        doubletime, sampletime, minvalue
    )
    return SimulationData(
        observedcases, 
        interventions, 
        Ns, 
        newexptdseedcases, 
        rng, 
        u0s, 
        betas, 
        mus, 
        deltas, 
        psis, 
        kappas, 
        id
    )
end

function Base.getproperty(obj::RenewalDiDDataUnlimitedPopn, k::Symbol)
    if k === :Ns 
        return nothing 
    else 
        return getfield(obj, k)
    end
end

function Base.show(io::IO, ::MIME"text/plain", d::AbstractRenewalDiDData)
    return _showabstractrenewaldiddata(io, d)
end

function _showabstractrenewaldiddata(io, d)
    _showtitle(io, d)
    print(io, "\n observedcases:  ")
    show(io, d.observedcases)
    print(io, "\n interventions:  ")
    show(io, d.interventions)
    print(io, "\n Ns:             ")
    _showns(io, d)
    print(io, "\n exptdseedcases: ")
    show(io, d.exptdseedcases)
    return nothing
end

function _showtitle(io, d::RenewalDiDData{S, T}) where {S, T}
    print(io, "RenewalDiDData{$S, $T}$(_printid(d))")
    return nothing
end

function _showtitle(io, d::RenewalDiDDataUnlimitedPopn{S, T}) where {S, T}
    print(io, "RenewalDiDDataUnlimitedPopn{$S, $T}$(_printid(d))")
    return nothing
end

function _showtitle(io, d::SimulationData{S, T}) where {S, T}
    print(io, "SimulationData{$S, $T}$(_printid(d))")
    return nothing
end

_printid(d) = d.id == "" ? "" : ", ($(d.id))"

_showns(io, d::AbstractRenewalDiDData) = show(io, d.Ns)
_showns(io, ::RenewalDiDDataUnlimitedPopn) = print(io, "unlimited")
