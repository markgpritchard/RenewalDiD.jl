# plot trace rank plot

function RenewalDiD.tracerankplot(args...; kwargs...)
    fig = Figure(; kwargs...)
    RenewalDiD.tracerankplot!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.tracerankplot!(
    fig::FigOrGridLayout, df; 
    ncols=nothing, nplots=nothing, kwargs...
)
    variables = _traceplotplotnames(df, ncols, nplots)
    return RenewalDiD.tracerankplot!(fig, df, variables; kwargs...)
end

function RenewalDiD.tracerankplot!(fig::FigOrGridLayout, df, variable; kwargs...)
    ax = Axis(fig[1, 1]; axiskws(; kwargs...)...)
    RenewalDiD.tracerankplot!(ax, df, variable; kwargs...)
    _traceplotylabels!(fig, df, variable; kwargs...)
    return nothing
end

function RenewalDiD.tracerankplot!(
    fig::FigOrGridLayout, df, variables::AbstractVector; 
    kwargs...
)
    axs = [Axis(fig[i, 1]; axiskws(; kwargs...)...) for i in eachindex(variables)]
    RenewalDiD.tracerankplot!(axs, df, variables; kwargs...)
    _traceplotylabels!(fig, df, variables; kwargs...)
    return nothing
end

function RenewalDiD.tracerankplot!(
    fig::FigOrGridLayout, df, variables::AbstractMatrix; 
    kwargs...
)
    axs = [
        ismissing(variables[i, j]) ? 
            nothing : 
            Axis(fig[i, (2 * j - 1)]; axiskws(; kwargs...)...)
        for i in axes(variables, 1), j in axes(variables, 2)
    ]
    RenewalDiD.tracerankplot!(axs, df, variables; kwargs...)
    for i in axes(variables, 1), j in axes(variables, 2) 
        _traceplotylabels!(fig, df, variables[i, j], i; col=(2 * (j-1)), kwargs...)
    end
    return nothing
end

function RenewalDiD.tracerankplot!(ax::Axis, df, variable; kwargs...)
    return _tracerankplot!(ax, df, variable; kwargs...)
end

function RenewalDiD.tracerankplot!(
    axs::AbstractArray{<:Union{<:Axis, Nothing}}, df, variables; 
    kwargs...
)
    return _tracerankplot!(axs, df, variables; kwargs...)
end

function _tracerankplot!(ax::Axis, df, variable::StringOrSymbol; binsize=10, kwargs...)
    rv = rankvalues(df, variable; binsize)
    for (i, chain) in enumerate(unique(df.chain))
        inds = findall(x -> x == chain, df.chain)
        lines!(ax, df.iteration[inds], rv[inds]; lineskws(; kwargs...)...)
    end
    return nothing
end

function _tracerankplot!(
    axs::AbstractArray{<:Union{<:Axis, Nothing}}, 
    df, 
    variables::AbstractArray{<:StringOrSymbolOrMissing}; 
    kwargs...
)
    for (i, ax) in enumerate(axs)
        _tracerankplot!(ax, df, variables[i]; kwargs...)
    end
    return nothing
end

function _tracerankplot!(
    axs, df, vs::X; 
    kwargs...
) where X <: Union{<:Integer, <:AbstractArray{<:Integer}}
    variables = _plotnames(df, vs)
    return _tracerankplot!(axs, df, variables; kwargs...)
end

_tracerankplot!(::Nothing, args...; kwargs...) = nothing
