# plot trace rank plot

function RenewalDiD.Plotting.tracerankplot(df, variable; binsize=10, kwargs...)
    fig = Figure(; kwargs...)
    ax = Axis(fig[1, 1]; axiskws(; kwargs...)...)
    tracerankplot!(ax, df, variable; binsize, kwargs...)
    return fig
end

function RenewalDiD.Plotting.tracerankplot!(ax, df, variable; binsize=10, kwargs...)
    _tracerankplot!(ax, df, variable; binsize=10, lineskws(kwargs...)...)
    return nothing
end

function _tracerankplot!(ax, df, variable; binsize=10, kwargs...)
    rv = rankvalues(df, variable; binsize)
    for (i, chain) in enumerate(unique(df.chain))
        inds = findall(x -> x == chain, df.chain)
        lines!(ax, df.iteration[inds], rv[inds]; kwargs...)
    end
    return nothing
end
