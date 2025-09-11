# The InterventionMatrix struct for storing times that interventions are implemented

# Structs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
    AbstractInterventionArray{T, N} <: AbstractArray{T, N}

Abstract type of array of intervention times.

Three abstract subtypes and three subtypes are defined:
- `InterventionVector{T} <: AbstractInterventionVector{T} <: AbstractInterventionArray{T, 1}`
- `InterventionMatrix{T} <: AbstractInterventionMatrix{T} <: AbstractInterventionArray{T, 2}`
- `InterventionArray{T} <: AbstractInterventionArray3{T} <: AbstractInterventionArray{T, 3}`
"""
abstract type AbstractInterventionArray{T, N} <: AbstractArray{T, N} end


"""
    AbstractInterventionVector{T} <: AbstractInterventionArray{T, 1}

Abstract type of vector of intervention times.

One type is defined: 
- `InterventionVector{T}`
"""
abstract type AbstractInterventionVector{T} <: AbstractInterventionArray{T, 1} end

"""
    InterventionVector{T}

Vector of intervention times.

# Fields
- `duration::Int`: duration of the analysis
- `offset::Int`: how much `rawstarttime` has been shifted (used for placebo interventions, 
    see examples)
- `rawstarttime::Union{<:Number, Nothing}`: time of intervention; `nothing`, and values 
    `≤1` or `>duration`, are treated equivalently as meaning that the intervention did not 
    occur during the analysis period
- `starttime::Int`: cleaned version of `rawstarttime` that is called when in use

# Constructor

    InterventionVector[{T}](duration::Number, rawstarttime::Union{<:Number, Nothing}; offset=0)

# Examples
```jldoctest
julia> InterventionVector(100, 25)
100-element InterventionVector{Int64}
 time │ 1 
──────┼───
    1 │ 0
    ⋮ │ ⋮
   25 │ 1
    ⋮ │ ⋮
  100 │ 1


julia> InterventionVector{Float64}(100, nothing)
100-element InterventionVector{Float64}
 time │   1 
──────┼─────
    1 │ 0.0
    ⋮ │   ⋮
  100 │ 0.0


julia> InterventionVector(100, 50; offset=-7)
100-element InterventionVector{Int64}
 time │ 1 
──────┼───
    1 │ 0
    ⋮ │ ⋮
   43 │ 1
    ⋮ │ ⋮
  100 │ 1


julia> InterventionVector(100, nothing) ==
       InterventionVector(100, 1) ==
       InterventionVector(100, 101) ==
       InterventionVector(100, 50; offset=51)
true
```
""" 
struct InterventionVector{T} <: AbstractInterventionVector{T}
    duration::Int
    offset::Int
    rawstarttime::Union{<:Number, Nothing}
    starttime::Int

    function InterventionVector{T}(
        duration::Number, rawstarttime::S;
        offset=0,
        mutewarnings=nothing,  # not used but allowed for consistency with InterventionMatrix
    ) where {S <: Union{<:Number, Nothing}, T <: Any}
        starttime = _interventionstarttime(rawstarttime, duration, offset)
        return new{T}(
            convert(Int, duration), convert(Int, offset), rawstarttime, starttime
        )
    end
end

InterventionVector(args...; kwargs...) = InterventionVector{Int}(args...; kwargs...)

"""
    InterventionVector(v::AbstractInterventionVector; offset)

Add an offset to an `InterventionVector`.

# Examples
```jldoctest
julia> v1 = InterventionVector(100, 50);

julia> InterventionVector(v1; offset=7)
100-element InterventionVector{Int64}
 time │ 1 
──────┼───
    1 │ 0
    ⋮ │ ⋮
   57 │ 1
    ⋮ │ ⋮
  100 │ 1
```
"""
function InterventionVector(v::AbstractInterventionVector; offset, kwargs...) 
    return _addoffsettointerventionarray(v, offset; kwargs...)
end

## InterventionMatrix

"""
    AbstractInterventionMatrix{T} <: AbstractInterventionArray{T, 2}

Abstract type of matrix of intervention times.

One type is defined: 
- `InterventionMatrix{T}`
"""
abstract type AbstractInterventionMatrix{T} <: AbstractInterventionArray{T, 2} end

"""
    InterventionMatrix{T}

Matrix of intervention times.

# Fields
- `duration::Int`: duration of the analysis
- `offset::Int`: how much `rawstarttime` has been shifted (used for placebo interventions, see 
    examples)
- `rawstarttimes::Vector{<:Union{<:Integer, Nothing}}`: times of interventions for each 
    group as entered by user; `nothing`, and values `≤1` or `>duration`, are treated 
    equivalently as meaning that the intervention did not occur during the analysis period
- `starttimes::Vector{Int}`: cleaned version of `rawstarttimes` that is called when in use

# Constructors

    InterventionMatrix[{T}](duration, rawstarttimes; mutewarnings=nothing, offset=0)
    InterventionMatrix[{T}](duration, rawstarttimes...; mutewarnings=nothing, offset=0)

`rawstarttimes` can be supplied as a vector of start times for each group or as separate 
    arguments.    

Warnings are provided if no groups have an intervention or all groups have an intervention
    before `duration`, unless `mutewarnings==true`.

# Examples
```jldoctest; filter=r"RenewalDiD*.*1045"
julia> InterventionMatrix(100, [25, 50, 200, 0, nothing])
100×5 InterventionMatrix{Int64}
 time │ 1  2  3  4  5 
──────┼───────────────
    1 │ 0  0  0  0  0
    ⋮ │ ⋮  ⋮  ⋮  ⋮  ⋮
   25 │ 1  0  0  0  0
    ⋮ │ ⋮  ⋮  ⋮  ⋮  ⋮
   50 │ 1  1  0  0  0
    ⋮ │ ⋮  ⋮  ⋮  ⋮  ⋮
  100 │ 1  1  0  0  0


julia> InterventionMatrix{Bool}(100, 0, 50, 95, nothing; offset=7)
100×4 InterventionMatrix{Bool}
 time │     1      2      3      4 
──────┼────────────────────────────
    1 │ false  false  false  false
    ⋮ │     ⋮      ⋮      ⋮      ⋮
    7 │  true  false  false  false
    ⋮ │     ⋮      ⋮      ⋮      ⋮
   57 │  true   true  false  false
    ⋮ │     ⋮      ⋮      ⋮      ⋮
  100 │  true   true  false  false


julia> InterventionMatrix(100, [25, 50, 100]);
┌ Warning: All groups in InterventionMatrix have intervention before end of duration
└ @ RenewalDiD yourpath.jl:1045
```
""" 
struct InterventionMatrix{T} <: AbstractInterventionMatrix{T}
    duration::Int
    offset::Int
    rawstarttimes::Vector{<:Union{<:Integer, Nothing}}
    starttimes::Vector{Int}

    function InterventionMatrix{T}(
        duration::Number, rawstarttimes::Vector{<:Union{<:Number, Nothing}};
        mutewarnings=nothing, offset=0,
    ) where T
        starttimes = _generateinterventionstarttimes(rawstarttimes, duration, offset)
        _interventionarrayconstructionwarnings(
            InterventionMatrix, duration, starttimes, mutewarnings
        )
        return new{T}(
            convert(Int, duration), convert(Int, offset), rawstarttimes, starttimes
        )
    end
end

InterventionMatrix(args...; kwargs...) = InterventionMatrix{Int}(args...; kwargs...)

function InterventionMatrix{T}(duration::Number, args...; kwargs...) where T
    return InterventionMatrix{T}(duration, [args...]; kwargs...)
end

function InterventionMatrix{T}(::Number, v::AbstractVector; kwargs...) where T
    # throw error and avoid a loop nesting vectors within vectors
    throw(_interventionarrayunusablevectorerror(v))
    return nothing
end

"""
    InterventionMatrix(M::AbstractInterventionMatrix; offset)

Add an offset to an `InterventionMatrix`.

# Examples
```jldoctest
julia> M = InterventionMatrix(100, [50, nothing]);

julia> InterventionMatrix(M; offset=7)
100×2 InterventionMatrix{Int64}
 time │ 1  2 
──────┼──────
    1 │ 0  0
    ⋮ │ ⋮  ⋮
   57 │ 1  0
    ⋮ │ ⋮  ⋮
  100 │ 1  0
```
"""
function InterventionMatrix(M::AbstractInterventionMatrix; offset, kwargs...)
    return _addoffsettointerventionarray(M, offset; kwargs...)
end

const InterventionVecOrMat{T} = Union{
    <:AbstractInterventionVector{T}, <:AbstractInterventionMatrix{T}
}

## Array of multiple interventions

"""
    AbstractInterventionArray3{T} <: AbstractInterventionArray{T, 3}

Abstract type of 3-dimensional array of intervention times.

One type is defined: 
- `InterventionArray{T}`
"""
abstract type AbstractInterventionArray3{T} <: AbstractInterventionArray{T, 3} end

"""
    InterventionArray{T}

Array of intervention times.

The third dimension represents different interventions. These can either be alternative
    interventions that were used, or offset times from the primary intervention.

# Fields
- `duration::Int`: duration of the analysis
- `offset::Vector{Int}`: how much the `rawstarttime` have been shifted for each type of 
    intervention (used for placebo interventions, see examples)
- `rawstarttimes::VecOrMat{<:Union{<:Integer, Nothing}}`: times of interventions for each 
    group as entered by user; `nothing`, and values `≤1` or `>duration`, are treated 
    equivalently as meaning that the intervention did not occur during the analysis period
- `starttimes::Matrix{Int}`: cleaned version of `rawstarttimes` that is called when in use

# Constructors

    InterventionArray[{T}](duration::Number, rawstarttimes::Matrix{<:Union{<:Number, Nothing}}; offset=0)
    InterventionArray[{T}](duration::Number, interventions::AbstractVector{<:AbstractVector}; offset=0)
    InterventionArray[{T}](duration::Number, intervention1::Vector, intervention2::Vector; offset=0)
    InterventionArray[{T}](duration::Number, intervention::Vector; offset=0)

Warnings are provided if no groups have an intervention or all groups have an intervention
    before the end of duration, unless the keyword argument `mutewarnings=true` is supplied.

If only a single `intervention` is supplied, a vector of `offset` values can be used to 
    generate an array of placebo intervention times.

# Examples
```jldoctest
julia> M = hcat([2, 3, nothing], [6, nothing, nothing]);

julia> InterventionArray(10, M)
"10×3×2 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    2 │ 1  0  0
    3 │ 1  1  0
    ⋮ │ ⋮  ⋮  ⋮
   10 │ 1  1  0

[:, :, 2] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
    6 │ 1  0  0
    ⋮ │ ⋮  ⋮  ⋮
   10 │ 1  0  0


julia> InterventionArray(10, M) ==
       InterventionArray(10, [[2, 3, nothing], [6, nothing, nothing]]) ==
       InterventionArray(10, [2, 3, nothing], [6, nothing, nothing])
true

julia> InterventionArray(10, M; offset=7)
"10×3×2 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
    9 │ 1  0  0
   10 │ 1  1  0

[:, :, 2] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
   10 │ 0  0  0


julia> InterventionArray(100, [2, 3, nothing]; offset=-7:7:7)
"100×3×3 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
  100 │ 0  0  0

[:, :, 2] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    2 │ 1  0  0
    3 │ 1  1  0
    ⋮ │ ⋮  ⋮  ⋮
  100 │ 1  1  0

[:, :, 3] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
    9 │ 1  0  0
   10 │ 1  1  0
    ⋮ │ ⋮  ⋮  ⋮
  100 │ 1  1  0
```
""" 
struct InterventionArray{T} <: AbstractInterventionArray3{T}
    duration::Int
    offset::Vector{Int}
    rawstarttimes::VecOrMat{<:Union{<:Integer, Nothing}}
    starttimes::Matrix{Int}

    function InterventionArray{T}(
        duration::Number, 
        offset::AbstractVector{<:Number}, 
        rawstarttimes::Matrix{<:Union{<:Number, Nothing}};
        mutewarnings=nothing,
    ) where T
        length(offset) == size(rawstarttimes, 2) || throw(ArgumentError("to do"))
        starttimes = _generateinterventionstarttimes(rawstarttimes, duration, offset)
        _interventionarrayconstructionwarnings(
            InterventionArray, duration, starttimes, mutewarnings
        )
        return new{T}(
            convert(Int, duration), 
            convert.(Int, offset), 
            rawstarttimes, 
            starttimes
        )
    end
end

"""
    InterventionArray(A::InterventionArray{T}; offset=0, kwargs...)
    InterventionArray(A::InterventionMatrix{T}; offset=[0], kwargs...)

Add an offset to an `InterventionArray` or a vector of offsets to an `InterventionMatrix`.

# Examples
```jldoctest
julia> A = InterventionArray(10, [2, 3, nothing], [6, nothing, nothing]);

julia> InterventionArray(A; offset=7)
"10×3×2 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
    9 │ 1  0  0
   10 │ 1  1  0

[:, :, 2] =
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
   10 │ 0  0  0

julia> M = InterventionMatrix(10, [2, nothing]);

julia> InterventionArray(M; offset=-7:7:7)
"10×2×3 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    ⋮ │ ⋮  ⋮
   10 │ 0  0

[:, :, 2] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    2 │ 1  0
    ⋮ │ ⋮  ⋮
   10 │ 1  0

[:, :, 3] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    ⋮ │ ⋮  ⋮
    9 │ 1  0
   10 │ 1  0
```
"""
function InterventionArray(A::InterventionArray{T}; offset=0, kwargs...) where T
    return __InterventionArray(T, A, offset; kwargs...)
end

function InterventionArray(A::InterventionMatrix{T}; offset=[0], kwargs...) where T
    return __InterventionArray(T, A, offset; kwargs...)
end

InterventionArray(args...; kwargs...) = InterventionArray{Int}(args...; kwargs...)

function InterventionArray{T}(args...; kwargs...) where T 
    return _InterventionArray(T, args...; kwargs...)
end

function _InterventionArray(
    T, duration::Number, intervention1::VecOrMat, intervention2::Vector; 
    kwargs...
) 
    return _InterventionArray(T, duration, hcat(intervention1, intervention2); kwargs...)
end

function _InterventionArray(
    T, duration::Number, intervention1::VecOrMat, intervention2::Vector, args...; 
    kwargs...
) 
    # loop through intervention vectors until a single matrix of intervention start times 
    return _InterventionArray(
        T, duration, hcat(intervention1, intervention2), args...; 
        kwargs...
    )
end

function _InterventionArray(
    T, duration::Number, v::AbstractVector{<:AbstractVector}; 
    kwargs...
) 
    # splat vector of vectors so the above versions can concatenate them into a matrix
    return _InterventionArray(T, duration, v...; kwargs...)
end

function _InterventionArray(
    T, duration::Number, interventions::AbstractVector; 
    offset=0, kwargs...
) 
    return __InterventionArray(T, duration, offset, interventions; kwargs...)
end

function _InterventionArray(
    T, duration::Number, interventions::AbstractMatrix; 
    offset=0, kwargs...
) 
    return __InterventionArray(T, duration, offset, interventions; kwargs...)
end

function __InterventionArray(
    T, duration::Number, offset::Number, interventions::AbstractMatrix; 
    kwargs...
) 
    # one offset and many sets of interventions -- make offsets match interventions
    v = [offset for _ in axes(interventions, 2)]
    return InterventionArray{T}(duration, v, interventions; kwargs...)
end
#=
function _InterventionArray(T, A::InterventionArray; offset=0, kwargs...) 
    return __InterventionArray(T, A, offset; kwargs...)
end
=#
function __InterventionArray(
    T, duration::Number, offset::Number, interventions::AbstractVector; 
    kwargs...
) 
    # only one offset and only one set of interventions -- essentially a 3-dimensional matrix
    return InterventionArray{T}(duration, [offset], [interventions;; ]; kwargs...)
end

function __InterventionArray(
    T, duration::Number, offset::AbstractVector, interventions::AbstractVector; 
    kwargs...
) 
    # multiple offsets and one set of interventions -- make interventions match offsets 
    M = hcat([interventions for _ in eachindex(offset)]...)
    return __InterventionArray(T, duration, offset, M; kwargs...)
end

function __InterventionArray(
    T, duration::Number, offset::AbstractVector, interventions::AbstractMatrix; 
    kwargs...
) 
    return InterventionArray{T}(duration, offset, interventions; kwargs...)
end

function __InterventionArray(T, A::InterventionArray, offset::Number; kwargs...) 
    v = [offset for _ in axes(A, 3)]
    return __InterventionArray(T, A, v; kwargs...)
end

function __InterventionArray(T, A::InterventionArray, offset::AbstractVector; kwargs...) 
    A isa InterventionArray{T} || throw(ArgumentError("to do"))
    return _addoffsettointerventionarray(A, offset; kwargs...)
end

function __InterventionArray(T, A::AbstractMatrix, offset::AbstractVector; kwargs...) 
    A isa AbstractMatrix{T} || throw(ArgumentError("to do"))
    return _addoffsettointerventionarray(A, offset; kwargs...)
end

## Functions to generate AbstractInterventionArray

function _generateinterventionstarttimes(rawstarttimes, duration, offset)
    starttimes = zeros(Int, size(rawstarttimes)...)
    _generateinterventionstarttimes!(starttimes, rawstarttimes, duration, offset)
    return starttimes
end

function _generateinterventionstarttimes!(
    starttimes::AbstractVector, rawstarttimes::AbstractVector, duration, offset::Number
)
    for (i, t) in enumerate(rawstarttimes)
        starttimes[i] = _interventionstarttime(t, duration, offset)
    end
    return nothing
end

function _generateinterventionstarttimes!(
    starttimes::Matrix, rawstarttimes::Matrix, duration, offset::AbstractVector
)
    for k in axes(rawstarttimes, 2)
        _generateinterventionstarttimes!(
            (@view starttimes[:, k]), (@view rawstarttimes[:, k]), duration, offset[k]
        )
    end
    return nothing
end

function _interventionstarttime(t::Number, duration, offset)
    if t + offset <= 1 || t + offset > duration
        return convert(Int, duration + 1)
    else
        return convert(Int, t + offset)
    end
end

_interventionstarttime(::Nothing, duration, ::Any) = convert(Int, duration + 1)

function _addoffsettointerventionarray(
    A::S, offset; 
    kwargs...
) where S <: AbstractInterventionArray
    totalnewoffset = A.offset .+ offset
    return S(A.duration, A.rawstarttimes; offset=totalnewoffset, kwargs...)
end

function _addoffsettointerventionarray(
    A::InterventionMatrix{T}, offset::AbstractVector; 
    kwargs...
) where T
    totalnewoffset = A.offset .+ offset
    rawstarttimesmatrix = repeat([A.rawstarttimes; ]; outer=(1, length(offset)))
    return InterventionArray{T}(
        A.duration, rawstarttimesmatrix; 
        offset=totalnewoffset, kwargs...
    )
end

## reproduce AbstractInterventionArray with different duration

"""
    modifiedduration(A, duration)

Reproduce an `AbstractInterventionArray` with different duration.

# Examples
```jldoctest
julia> v = InterventionVector(10, 7);

julia> modifiedduration(v, 5)
5-element InterventionVector{Int64}
 time │ 1 
──────┼───
    1 │ 0
    ⋮ │ ⋮
    5 │ 0


julia> M = InterventionMatrix(10, [2, nothing]);

julia> modifiedduration(M, 5)
5×2 InterventionMatrix{Int64}
 time │ 1  2 
──────┼──────
    1 │ 0  0
    2 │ 1  0
    ⋮ │ ⋮  ⋮
    5 │ 1  0

    
julia> A = InterventionArray(10, [2, nothing], [6, nothing]);

julia> modifiedduration(A, 5)
"5×2×2 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    2 │ 1  0
    ⋮ │ ⋮  ⋮
    5 │ 1  0

[:, :, 2] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    ⋮ │ ⋮  ⋮
    5 │ 0  0
```
"""
function modifiedduration(A::T, duration; kwargs...) where T <: AbstractInterventionArray 
    offset = A.offset 
    rawstarttimes = A.rawstarttimes
    return T(duration, rawstarttimes; offset)
end


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Properties of `InterventionMatrix`

function Base.getproperty(v::AbstractInterventionVector, k::Symbol)
    if k === :starttimes
        return getfield(v, :starttime)
    elseif k === :rawstarttimes
        return getfield(v, :rawstarttime)
    else 
        return getfield(v, k)
    end
end

Base.size(v::AbstractInterventionVector) = (v.duration, )
Base.size(M::AbstractInterventionMatrix) = (M.duration, length(M.starttimes))
Base.size(A::AbstractInterventionArray3) = (A.duration, size(A.starttimes)...)

# `Base.getindex` defined for Cartesian Integer indexes; for all other indexing, functions 
# in `Base` are used until a tuple of integers are passed 
function Base.getindex(v::AbstractInterventionVector{T}, t::Integer) where T
    0 < t <= _duration(v) || throw(BoundsError(v, t))
    return T(t >= v.starttime)
end

function Base.getindex(M::AbstractInterventionMatrix{T}, t::Integer, x::Integer) where T
    0 < t <= _duration(M) || throw(BoundsError(M, [t, x]))
    0 < x <= length(M.starttimes) || throw(BoundsError(M, [t, x]))
    return T(t >= M.starttimes[x])
end

function Base.getindex(
    A::AbstractInterventionArray3{T}, t::Integer, x::Integer, z::Integer
) where T
    0 < t <= _duration(A) || throw(BoundsError(A, [t, x, z]))
    0 < x <= size(A.starttimes, 1) || throw(BoundsError(A, [t, x, z]))
    0 < z <= size(A.starttimes, 2) || throw(BoundsError(A, [t, x, z]))
    return T(t >= A.starttimes[x, z])
end

_duration(A::AbstractInterventionArray) = A.duration
_offset(A::AbstractInterventionArray) = A.offset
_starttimes(A::AbstractInterventionVector) = A.starttime
_starttimes(A::AbstractInterventionArray) = A.starttimes

## Intervention times for plot 

function _interventionstarttimes(M::AbstractMatrix)
    return [_interventionstarttimes(M, i) for i in axes(M, 2)]
end

function _interventionstarttimes(M::AbstractInterventionMatrix, i)
    return M.starttimes[i] > _duration(M) ? nothing : M.starttimes[i]
end  # another version of this function is in `processparameters.jl`

## Concatenate 

function Base.cat(A::AbstractInterventionArray, args...; dims, kwargs...)
    return _cat(dims, A, args...; kwargs...)
end

_cat(dims::Number, args...; kwargs...) = _cat(Val(dims), args...; kwargs...)

function _cat(::Val{1}, args...; kwargs...)
    throw(_catd1error())
end

function _cat(::Val{2}, A::InterventionVecOrMat{T}, args...; mutewarnings=nothing) where T
    d1 = A.duration
    o1 = A.offset
    for a in args
        a isa InterventionVecOrMat{T} || throw(_cattypeerror("InterventionMatrix{$T}", a, 2))
        a.duration == d1 || throw(_durationmismatcherror(d1, a))
        a.offset == o1 || throw(_unequaloffseterror(o1, a))
    end
    rawstarttimes = vcat(A.rawstarttimes, [a.rawstarttimes for a in args]...)
    return InterventionMatrix{T}(d1, rawstarttimes; offset=o1, mutewarnings)
end

function _cat(::Val{2}, A::AbstractInterventionArray3{T}, args...) where T
    d1 = A.duration
    for a in args
        a isa AbstractInterventionArray3{T} || 
            throw(_cattypeerror("AbstractInterventionArray3{$T}", a, 2))
        a.duration == d1 || throw(_durationmismatcherror(d1, a))
    end
    rawstarttimes = vcat(A.rawstarttimes, [a.rawstarttimes for a in args]...)
    return InterventionArray{T}(d1, rawstarttimes; )
end

# combined type only used in the following function signature
const _MatOrA3{T} = Union{<:AbstractInterventionMatrix{T}, <:AbstractInterventionArray3{T}}

function _cat(::Val{3}, A::_MatOrA3{T}, args...; kwargs...) where T
    d1 = A.duration
    for a in args
        a isa _MatOrA3{T} || throw(_cattypeerror("InterventionArray{$T}", a, 3))
        a.duration == d1 || throw(_durationmismatcherror(d1, a))
        size(a, 2) == size(A, 2) || throw(_groupnumbermismatch(A, a))
    end
    offset = vcat(A.offset, [a.offset for a in args]...)
    rawstarttimes = hcat(A.rawstarttimes, [a.rawstarttimes for a in args]...)
    return InterventionArray{T}(d1, offset, rawstarttimes; kwargs...)
end

"""
    interventioncat(args...)

Concatenate `AbstractInterventionMatrix` or `AbstractInterventionArray3` to give multiple 
    interventions.

Can also use `Base.cat(args...; dims=3)`.

# Examples
```jldoctest
julia> M1 = InterventionMatrix(10, [2, nothing]);

julia> M2 = InterventionMatrix(10, [2, 6]; offset=5);

julia> interventioncat(M1, M2)
"10×2×2 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    2 │ 1  0
    ⋮ │ ⋮  ⋮
   10 │ 1  0

[:, :, 2] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    ⋮ │ ⋮  ⋮
    7 │ 1  0
    ⋮ │ ⋮  ⋮
   10 │ 1  0


julia> A = InterventionArray(10, [2, nothing], [6, nothing]);

julia> interventioncat(A, M1, M2)
"10×2×4 InterventionArray{Int64}"
[:, :, 1] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    2 │ 1  0
    ⋮ │ ⋮  ⋮
   10 │ 1  0

[:, :, 2] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    ⋮ │ ⋮  ⋮
    6 │ 1  0
    ⋮ │ ⋮  ⋮
   10 │ 1  0

[:, :, 3] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    2 │ 1  0
    ⋮ │ ⋮  ⋮
   10 │ 1  0

[:, :, 4] =
 time │ 1  2 
──────┼──────
    1 │ 0  0
    ⋮ │ ⋮  ⋮
    7 │ 1  0
    ⋮ │ ⋮  ⋮
   10 │ 1  0
```
"""
interventioncat(args...; kwargs...) = _cat(Val{3}(), args...; kwargs...) 

## Show

function Base.show(io::IO, ::MIME"text/plain", A::InterventionVecOrMat) 
    combinedstrings = _showcombinedstrings(A, nothing)
    @static if pkgversion(PrettyTables).major == 2
        return pretty_table(
            io, combinedstrings; 
            header=_showheader(A),
            hlines=[1], 
            show_row_number=false, 
            title=summary(A), 
            vlines=[1], 
        )
    else 
        return pretty_table(
            io, combinedstrings; 
            column_labels=_showheader(A),
            row_labels=nothing, 
            title=summary(A), 
            table_format=TextTableFormat( ;
                @text__no_horizontal_lines,
                horizontal_line_after_column_labels=true,
                vertical_line_at_beginning=false,
                vertical_lines_at_data_columns=[1], 
                vertical_line_after_data_columns=false,
            ),
        )
    end
end

function Base.show(io::IO, ::MIME"text/plain", A::AbstractInterventionArray3) 
    show(io, summary(A))
    if size(A, 3) <= 6
        for k in 1:6
            _showoneintervention(io, A, k)
        end
    else
        for k in 1:3
            _showoneintervention(io, A, k)
        end
        print("\n;;; …\n")
        s = size(A, 3)
        for k in s-2:s
            _showoneintervention(io, A, k)
        end
    end
    return nothing
end

# short version used when `InterventionMatrix` is in an `AbstractRenewalDiDData` 
function Base.show(io::IO, M::InterventionVecOrMat) 
    show(io, collect(M))
    print(io, " {duration $(M.duration), starttimes [")
    join(io, _showliststarttimes(M), ", ")
    print(io, "]}")
    return nothing
end

Base.show(io::IO, A::AbstractInterventionArray3) = show(io, collect(collect(A)))

function _showoneintervention(io, A, k)
    k > size(A, 3) && return nothing 
    print(io, "\n[:, :, $k] =\n")
    combinedstrings = _showcombinedstrings(A, k)
    @static if pkgversion(PrettyTables).major == 2
        return pretty_table(
            io, combinedstrings; 
            header=_showheader(A),
            hlines=[1], 
            show_row_number=false, 
            vlines=[1], 
        )
    else 
        return pretty_table(
            io, combinedstrings; 
            column_labels=_showheader(A),
            row_labels=nothing, 
            table_format=TextTableFormat( ;
                @text__no_horizontal_lines,
                horizontal_line_after_column_labels=true,
                vertical_line_at_beginning=false,
                vertical_lines_at_data_columns=[1], 
                vertical_line_after_data_columns=false,
            ),
        )
    end

end

_showliststarttimes(M) = [x > M.duration ? nothing : x for x in M.starttimes]

function _showtimes(M, k)
    uniquestarttimes = __uniquestarttimes(M, k)
    unsortedshowtimes = unique([uniquestarttimes; 1; _duration(M)])
    showtimes = sort(unsortedshowtimes)
    filter!(x -> x <= _duration(M), showtimes)
    return showtimes
end

__uniquestarttimes(M::InterventionVecOrMat, ::Nothing) = unique(M.starttimes)
__uniquestarttimes(M::AbstractInterventionArray3, k) = unique(@view M.starttimes[:, k])

function _showcontents(M::InterventionVecOrMat{T}, showtimes, ::Nothing) where T
    C = zeros(T, length(showtimes), size(M, 2))
    for (i, t) in enumerate(showtimes), g in axes(M, 2)
        C[i, g] = M[t, g]
    end 
    return C 
end

function _showcontents(A::AbstractInterventionArray3{T}, showtimes, k) where T 
    C = zeros(T, length(showtimes), size(A, 2))
    for (i, t) in enumerate(showtimes), g in axes(A, 2)
        C[i, g] = A[t, g, k]
    end 
    return C 
end

function _showstrings(M::AbstractInterventionArray{T}, k) where T
    showtimes = _showtimes(M, k)
    contents = _showcontents(M, showtimes, k)
    stringtimes = ["1"]
    stringcontents = _showstringmatrixrepeatedrow(M, "$(zero(T))")

    for (i, t) in enumerate(showtimes) 
        t == 1 && continue 

        if t > showtimes[i-1] + 1
            push!(stringtimes, "⋮")
            stringcontents = vcat(stringcontents, _matrixvdots(stringcontents))
        end

        push!(stringtimes, "$t")
        r = _showstringmatrixrow(M)
        for g in axes(r, 2)
            r[1, g] = "$(contents[i, g])"
        end
        stringcontents = vcat(stringcontents, r)
    end

    return (stringtimes, stringcontents)
end

_showstringmatrixrow(M) = Matrix{String}(undef, 1, size(M, 2))

function _showstringmatrixrepeatedrow(M, x)
    r = _showstringmatrixrow(M)
    for g in axes(r, 2)
        r[1, g] = x
    end
    return r
end

_matrixvdots(M) = _showstringmatrixrepeatedrow(M, "⋮")

function _showcombinedstrings(M, k)
    stringtimes, stringcontents = _showstrings(M, k)
    return hcat(stringtimes, stringcontents)
end

_showheader(M) = ["time"; ["$x" for x in axes(M, 2)]]


# Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _interventionarrayconstructionwarnings(T, duration, starttimes, mutewarnings::Bool)
    if mutewarnings 
        return nothing 
    else 
        return _interventionarrayconstructionwarnings(T, duration, starttimes, nothing)
    end
end

function _interventionarrayconstructionwarnings(T, duration, starttimes, ::Nothing)
    minimum(starttimes) > duration && _nointerventionwarning(T)
    maximum(starttimes) <= duration && _allinterventionwarning(T)
    return nothing
end

_nointerventionwarning(T) = @warn "$T with no intervention in any group"

function _allinterventionwarning(T)
    @warn "All groups in $T have intervention before end of duration"
    return nothing
end


# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _catd1error()
    m = "concatenating `AbstractInterventionArray`s on dimension 1 is intentionally \
        unsupported"
    return ArgumentError(m)
end

function _cattypeerror(T, a, dim)
    m = "cannot concatenate a $(typeof(a)) to an $T on dimension $dim"
    return ArgumentError(m)
end

function _durationmismatcherror(d1, a)
    return DimensionMismatch("mismatch in duration (expected $d1 got $(a.duration))")
end

function _groupnumbermismatch(A, a)
    m = "mismatch in number of groups (expected $(size(A, 2)) got $(size(a, 2)))"
    DimensionMismatch(m)
end

function _interventionarrayunusablevectorerror(v)
    return ArgumentError("$v, vector received that could not be used as start times")
end

function _unequaloffseterror(o1, a)
    m = "when concatenating intevention arrays on dimension 2, all must have equal offset \
        (expected $o1, got $(a.offset))"
    return ArgumentError(m)
end
