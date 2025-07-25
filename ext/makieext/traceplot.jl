# plot trace rank plot

function RenewalDiD.Plotting.traceplot(df, variable; kwargs...)
    fig = Figure(; kwargs...)
    ax = Axis(fig[1, 1]; axiskws(; kwargs...)...)
    traceplot!(ax, df, variable; kwargs...)
    return fig
end

function RenewalDiD.Plotting.traceplot!(ax, df, variable; kwargs...)
    _traceplot!(ax, df, variable; lineskws(kwargs...)...)
    return nothing
end

function _traceplot!(ax, df, variable; kwargs...)
    for (i, chain) in enumerate(unique(df.chain))
        inds = findall(x -> x == chain, df.chain)
        lines!(ax, df.iteration[inds], getproperty(df, variable)[inds]; kwargs...)
    end
    return nothing
end
