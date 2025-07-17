# The InterventionMatrix struct for storing times that interventions are implemented

## Structs and constructors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abstract type AbstractInterventionsArray{T, N} <: AbstractArray{T, N} end

"""
    InterventionMatrix{T}

Matrix of intervention times.

It is expected that at least one group is treated and at least one group remains untreated 
at the end of the analysis. Untreated groups can be represented by a "start time" that is 
negative, greater than `duration`, or `nothing`. 

`T` is optional and defaults to `Int`.

# Fields
- `duration::Int`: duration of the analysis
- `rawstarttimes::Vector{<:Union{<:Integer, Nothing}}`: times of interventions for each 
    group (these are the values provided by the user when constructing the matrix; they are 
    retained for situations in which we wish to manipulate the start times)
- `starttimes::Vector{Int}`: times of interventions for each group (these are cleaned to a 
    vector of `Int` for efficiency)

# Constructors
Constructors accept only `duration` and `rawstarttimes`. `rawstarttimes` can be provided as 
    a vector of values or as separate arguments. 


# Examples
```
julia> M = InterventionMatrix(10, 2, 3, nothing)
10×3 InterventionMatrix{Int64}
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    2 │ 1  0  0
    3 │ 1  1  0
    ⋮ │ ⋮  ⋮  ⋮
   10 │ 1  1  0


julia> M == InterventionMatrix(10, [2, 3, nothing]) == InterventionMatrix(10, [2, 3, 100])
true

julia> M[5, 1]  # rows without new interventions are hidden but are still accessible
1
```
"""
struct InterventionMatrix{T} <: AbstractInterventionsArray{T, 2}
    duration::Int
    rawstarttimes::Vector{<:Union{<:Integer, Nothing}}
    starttimes::Vector{Int}

    function InterventionMatrix{T}(
        duration::Int, rawstarttimes::Vector{<:Union{<:Integer, Nothing}};
        mutewarnings=nothing,
    ) where T
        starttimes = _interventionstarttimes(rawstarttimes, duration)
        minimum(starttimes) > duration && _nointerventionwarning(mutewarnings)
        maximum(starttimes) <= duration && _allinterventionwarning(mutewarnings)
        return new{T}(duration, rawstarttimes, starttimes)
    end
end

### Constructors

InterventionMatrix(args...; kwargs...) = InterventionMatrix{Int}(args...; kwargs...)

function InterventionMatrix{T}(duration, args...; kwargs...) where T
    return _interventionmatrix(T, duration, args...; kwargs...)
end

function _interventionmatrix(T::DataType, duration::Integer, args...; kwargs...)
    return _interventionmatrix(T, convert(Int, duration), args...; kwargs...)
end

function _interventionmatrix(T::DataType, duration::Int, args...; kwargs...)
    return _interventionmatrix(T, duration, [args...]; kwargs...)
end

function _interventionmatrix(T::DataType, duration::Int, rawstarttimes::Vector; kwargs...)
    return __interventionmatrix(T, duration, rawstarttimes; kwargs...)
end

function __interventionmatrix(
    T::DataType, duration::Int, rawstarttimes::Vector{<:Union{<:Integer, Nothing}}; 
    kwargs...
)
    return InterventionMatrix{T}(duration, rawstarttimes; kwargs...)
end

function _interventionstarttimes(rawstarttimes, duration)
    starttimes = zeros(Int, length(rawstarttimes))
    for (i, t) in enumerate(rawstarttimes)
        if isnothing(t) || t <= 1 || t > duration
            starttimes[i] = duration + 1
        else
            starttimes[i] = convert(Int, t)
        end
    end
    return starttimes
end


## Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Base.size(M::InterventionMatrix) = (M.duration, length(M.starttimes))

function Base.getindex(M::InterventionMatrix{T}, t, x) where T
    t <= _duration(M) || throw(BoundsError(M, [t, x]))
    x <= length(M.starttimes) || throw(BoundsError(M, [t, x]))
    return T(t >= M.starttimes[x])
end

_duration(M::InterventionMatrix) = M.duration

function _showtimes(M::InterventionMatrix)
    uniquestarttimes = unique(M.starttimes)
    unsortedshowtimes = unique([uniquestarttimes; 1; _duration(M)])
    showtimes = sort(unsortedshowtimes)
    filter!(x -> x <= _duration(M), showtimes)
    return showtimes
end

function _showcontents(M::InterventionMatrix{T}, showtimes) where T 
    C = zeros(T, length(showtimes), size(M, 2))
    for (i, t) in enumerate(showtimes), g in axes(M, 2)
        C[i, g] = M[t, g]
    end 
    return C 
end

function _showstrings(M::InterventionMatrix{T}) where T
    showtimes = _showtimes(M)
    contents = _showcontents(M, showtimes)
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

function _showcombinedstrings(M::InterventionMatrix)
    stringtimes, stringcontents = _showstrings(M)
    return hcat(stringtimes, stringcontents)
end

_showheader(M::InterventionMatrix) = ["time"; ["$x" for x in axes(M, 2)]]

function Base.show(io::IO, ::MIME"text/plain", M::InterventionMatrix) 
    combinedstrings = _showcombinedstrings(M)
    return pretty_table(
        io, combinedstrings; 
        header=_showheader(M),
        hlines=[1], 
        show_row_number=false, 
        title=summary(M), 
        vlines=[1], 
    )
end


## Warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function _nointerventionwarning(::Nothing)
    @warn "InterventionMatrix with no intervention in any group"
    return nothing
end

function _nointerventionwarning(mutewarnings::Bool)
    if mutewarnings 
        return nothing 
    else 
        return _nointerventionwarning(nothing)
    end
end

function _allinterventionwarning(::Nothing)
    @warn "All groups in InterventionMatrix have intervention before end of duration"
    return nothing
end

function _allinterventionwarning(mutewarnings::Bool)
    if mutewarnings 
        return nothing 
    else 
        return _allinterventionwarning(nothing)
    end
end
