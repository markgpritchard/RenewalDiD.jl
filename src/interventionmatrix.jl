# The InterventionMatrix struct for storing times that interventions are implemented

# Structs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abstract type AbstractInterventionsArray{T, N} <: AbstractArray{T, N} end

"""
    InterventionMatrix{T}

Matrix of intervention times.

# Fields
- `duration::Int`: duration of the analysis
- `rawstarttimes::Vector{<:Union{<:Integer, Nothing}}`: times of interventions for each 
    group as entered by user
- `starttimes::Vector{Int}`: cleaned version of `rawstarttimes` that is called when in use

# Constructors

    InterventionMatrix[{T}](duration, rawstarttimes; mutewarnings=nothing)
    InterventionMatrix[{T}](duration, rawstarttimes...; mutewarnings=nothing)

## Arguments

* `T`: type of the array output; optional, defaults to `Int`.
* `duration`: must be able to convert into an `Int` 
* `rawstarttimes`: a vector or separate arguments; each must be `nothing` or a number that
    can convert into `Int`. Starttimes that are entered as `nothing`, a value `≤ 1` or 
    `> duration`, are treated equivalently as meaning that the intervention did not occur 
    during the study period.

Warnings are provided if no groups have an intervention or all groups have an intervention
    before the end of duration, unless `mutewarnings==true`.

# Examples
```jldoctest
julia> InterventionMatrix(100, [25, 50, 200])
100×3 InterventionMatrix{Int64}
 time │ 1  2  3 
──────┼─────────
    1 │ 0  0  0
    ⋮ │ ⋮  ⋮  ⋮
   25 │ 1  0  0
    ⋮ │ ⋮  ⋮  ⋮
   50 │ 1  1  0
    ⋮ │ ⋮  ⋮  ⋮
  100 │ 1  1  0


julia> InterventionMatrix{Bool}(100, 50, 75, -1, nothing)
100×4 InterventionMatrix{Bool}
 time │     1      2      3      4 
──────┼────────────────────────────
    1 │ false  false  false  false
    ⋮ │     ⋮      ⋮      ⋮      ⋮
   50 │  true  false  false  false
    ⋮ │     ⋮      ⋮      ⋮      ⋮
   75 │  true   true  false  false
    ⋮ │     ⋮      ⋮      ⋮      ⋮
  100 │  true   true  false  false


julia> InterventionMatrix(100, [25, 50, 100]);
┌ Warning: All groups in InterventionMatrix have intervention before end of duration
└ @ RenewalDiD 
```
""" 
struct InterventionMatrix{T} <: AbstractInterventionsArray{T, 2}
    duration::Int
    rawstarttimes::Vector{<:Union{<:Integer, Nothing}}
    starttimes::Vector{Int}

    function InterventionMatrix{T}(
        duration::Number, rawstarttimes::Vector{<:Union{<:Number, Nothing}};
        mutewarnings=nothing,
    ) where T
        starttimes = _generateinterventionstarttimes(rawstarttimes, duration)
        minimum(starttimes) > duration && _nointerventionwarning(mutewarnings)
        maximum(starttimes) <= duration && _allinterventionwarning(mutewarnings)
        return new{T}(convert(Int, duration), rawstarttimes, starttimes)
    end
end

## Constructors

InterventionMatrix(args...; kwargs...) = InterventionMatrix{Int}(args...; kwargs...)

function InterventionMatrix{T}(duration::Number, args...; kwargs...) where T
    return InterventionMatrix{T}(duration, [args...]; kwargs...)
end

function _generateinterventionstarttimes(rawstarttimes, duration)
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


# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Properties of `InterventionMatrix`

Base.size(M::InterventionMatrix) = (M.duration, length(M.starttimes))

function Base.getindex(M::InterventionMatrix{T}, t, x) where T
    t <= _duration(M) || throw(BoundsError(M, [t, x]))
    x <= length(M.starttimes) || throw(BoundsError(M, [t, x]))
    return T(t >= M.starttimes[x])
end

_duration(M::AbstractInterventionsArray) = M.duration

## Intervention times for plot 

function _interventionstarttimes(M::AbstractMatrix)
    return [_interventionstarttimes(M, i) for i in axes(M, 2)]
end

function _interventionstarttimes(M::InterventionMatrix, i)
    return M.starttimes[i] > _duration(M) ? nothing : M.starttimes[i]
end  # another version of this function is in `processparameters.jl`

## Show

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

# short version used when `InterventionMatrix` is in an `AbstractRenewalDiDData` 
function Base.show(io::IO, M::InterventionMatrix) 
    show(io, collect(M))
    print(io, " {duration $(M.duration), starttimes [")
    join(io, _showliststarttimes(M), ", ")
    print(io, "]}")
    return nothing
end

function _showliststarttimes(M)
    return [x > M.duration ? nothing : x for x in M.starttimes]
end

function _showtimes(M)
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

function _showcombinedstrings(M)
    stringtimes, stringcontents = _showstrings(M)
    return hcat(stringtimes, stringcontents)
end

_showheader(M) = ["time"; ["$x" for x in axes(M, 2)]]


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
