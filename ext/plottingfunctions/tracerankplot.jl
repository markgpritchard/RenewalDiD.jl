# plot trace rank plot

function RenewalDiD.Plotting.tracerankplot(df, variable; binsize=10)
    fig = Figure()
    ax = Axis(fig[1, 1])
    tracerankplot!(ax, df, variable; binsize)
    return fig
end

function RenewalDiD.Plotting.tracerankplot!(ax, df, variable; binsize=10)
    rv = rankvalues(df, variable; binsize)
    for (i, chain) in enumerate(unique(df.chain))
        inds = findall(x -> x == chain, df.chain)
        lines!(ax, df.iteration[inds], rv[inds])
    end
    return nothing
end
