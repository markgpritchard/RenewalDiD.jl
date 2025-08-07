# plot trace rank plot

function RenewalDiD.traceplot(args...; kwargs...)
    fig = Figure(; kwargs...)
    RenewalDiD.traceplot!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.traceplot!(
    fig::FigOrGridLayout, df; 
    ncols=nothing, nplots=nothing, kwargs...
)
    variables = _traceplotplotnames(df, ncols, nplots)
    return RenewalDiD.traceplot!(fig, df, variables; kwargs...)
end

function RenewalDiD.traceplot!(fig::FigOrGridLayout, df, variable; kwargs...)
    ax = Axis(fig[1, 1]; axiskws(; kwargs...)...)
    RenewalDiD.traceplot!(ax, df, variable; kwargs...)
    _traceplotylabels!(fig, df, variable; kwargs...)
    return nothing
end

function RenewalDiD.traceplot!(
    fig::FigOrGridLayout, df, variables::AbstractVector; 
    kwargs...
)
    axs = [Axis(fig[i, 1]; axiskws(; kwargs...)...) for i in eachindex(variables)]
    RenewalDiD.traceplot!(axs, df, variables; kwargs...)
    _traceplotylabels!(fig, df, variables; kwargs...)
    return nothing
end

function RenewalDiD.traceplot!(
    fig::FigOrGridLayout, df, variables::AbstractMatrix; 
    kwargs...
)
    axs = [
        ismissing(variables[i, j]) ? 
            nothing : 
            Axis(fig[i, (2 * j - 1)]; axiskws(; kwargs...)...)
        for i in axes(variables, 1), j in axes(variables, 2)
    ]
    RenewalDiD.traceplot!(axs, df, variables; kwargs...)
    for i in axes(variables, 1), j in axes(variables, 2) 
        _traceplotylabels!(fig, df, variables[i, j], i; col=(2 * (j-1)), kwargs...)
    end
    return nothing
end

function RenewalDiD.traceplot!(ax::Axis, df, variable; kwargs...)
    return _traceplot!(ax, df, variable; kwargs...)
end

function RenewalDiD.traceplot!(
    axs::AbstractArray{<:Union{<:Axis, Nothing}}, df, variables; 
    kwargs...
)
    return _traceplot!(axs, df, variables; kwargs...)
end

_traceplotplotnames(df, ::Nothing, ::Nothing) = _plotnames(df)

function _traceplotplotnames(df, ::Nothing, nplots::Integer)
    pn = _plotnames(df)
    if length(pn) <= nplots 
        return pn 
    else 
        return pn[1:nplots]
    end
end

function _traceplotplotnames(df, ncols::Integer, nplots::Any)
    pn = _traceplotplotnames(df, nothing, nplots)
    return _traceplotplotnamesmatrix(df, ncols, pn)
end

function _traceplotplotnamesmatrix(df, ncols, pn)
    nrows = round(Int, length(pn) / ncols, RoundUp)
    variables = Matrix{Union{String, Missing}}(missing, nrows, ncols)
    for i in eachindex(pn)
        variables[i] = pn[i]
    end
    return variables 
end

function _traceplot!(ax::Axis, df, variable::StringOrSymbol; kwargs...)
    for (i, chain) in enumerate(unique(df.chain))
        inds = findall(x -> x == chain, df.chain)
        lines!(
            ax, df.iteration[inds], getproperty(df, variable)[inds]; 
            lineskws(; kwargs...)...
        )
    end
    return nothing
end

_traceplot!(::Axis, ::Any, ::Missing; kwargs...) = nothing

function _traceplot!(
    axs::AbstractArray{<:Union{<:Axis, Nothing}}, 
    df, 
    variables::AbstractArray{<:StringOrSymbolOrMissing}; 
    kwargs...
)
    for (i, ax) in enumerate(axs)
        _traceplot!(ax, df, variables[i]; kwargs...)
    end
    return nothing
end

function _traceplot!(
    axs, df, vs::X; 
    kwargs...
) where X <: Union{<:Integer, <:AbstractArray{<:Integer}}
    variables = _plotnames(df, vs)
    return _traceplot!(axs, df, variables; kwargs...)
end

_traceplot!(::Nothing, args...; kwargs...) = nothing

function _traceplotylabels!(
    fig::FigOrGridLayout, ::Any, vs::AbstractVector{<:StringOrSymbolOrMissing}; 
    kwargs...
)
    for (i, v) in enumerate(vs)
        _traceplotylabels!(fig, nothing, v, i; kwargs...)
    end
    return nothing
end

function _traceplotylabels!(fig::FigOrGridLayout, ::Any, v::Symbol, row=1; kwargs...)
    return _traceplotylabels!(fig, nothing, String(v), row; kwargs...)
end

function _traceplotylabels!(
    fig::FigOrGridLayout, ::Any, v::AbstractString, row=1; 
    col=0, fontsize=11.84, rotation=π/2, tellheight=false, kwargs...
)
    Label(fig[row, col], v; labelkws(; fontsize, rotation, tellheight, kwargs...)...)
    return nothing
end

function _traceplotylabels!(
    fig::FigOrGridLayout, df, vs::AbstractArray{<:Integer}; 
    kwargs...
) 
    variables = _plotnames(df, vs)
    return _traceplotylabels!(fig, df, variables; kwargs...)
end

function _traceplotylabels!(fig::FigOrGridLayout, df, v::Integer, row=1; kwargs...)
    variable = _plotnames(df, v)
    return _traceplotylabels!(fig, df, variable, row; kwargs...)
end

_traceplotylabels!(::FigOrGridLayout, ::Any, ::Missing, args...; kwargs...) = nothing
