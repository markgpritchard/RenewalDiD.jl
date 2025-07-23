# plot trace rank plot

function RenewalDiD.traceplot(df, variable)
    fig = Figure()
    ax = Axis(fig[1, 1])
    traceplot!(ax, df, variable)
    return fig
end

function RenewalDiD.traceplot!(ax, df, variable)
    for (i, chain) in enumerate(unique(df.chain))
        inds = findall(x -> x == chain, df.chain)
        lines!(ax, df.iteration[inds], getproperty(df, variable)[inds])
    end
    return nothing
end
