# plot output from model

function RenewalDiD.Plotting.plotmodeloutput(args...; kwargs...)
    fig = Figure(; kwargs...)
    plotmodeloutput!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.Plotting.plotmodeloutput!(
    gl::FigOrGridLayout, A::AbstractArray, t=automatic; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, A; axiskws(; kwargs...)...)
    _plotmodeloutput!(axs, A, t; kwargs...)
    linkaxes!(axs...)
    return nothing
end

function RenewalDiD.Plotting.plotmodeloutput!(
    axs::AbstractVector{Axis}, A::AbstractArray, t=automatic; 
    kwargs...
) 
    length(axs) == size(A, 2) || throw(ArgumentError("to do"))
    _plotmodeloutput!(axs, A, t; kwargs...)
    return nothing
end

_plotmodeloutputaxs(gl, A; kwargs...) = [Axis(gl[1, i]; kwargs...) for i in axes(A, 2)]

function _plotmodeloutput!(axs, A, ::Automatic; kwargs...)
    t = axes(A, 1)
    _plotmodeloutput!(axs, A, t; kwargs...)
    return nothing
end

function _plotmodeloutput!(
    axs, A::AbstractMatrix, t::AbstractVector; 
    color=Cycled(1), kwargs...
) 
    _plotmodeloutput!(axs, A, t, 1; color, kwargs...)
    return nothing
end

function _plotmodeloutput!(
    axs, A::AbstractArray{T, 3}, t::AbstractVector; 
    color=Cycled(1), kwargs...
) where T
    nquantiles = size(A, 3)
    isodd(size(A, 3)) || throw(ArgumentError("to do"))
    _plotmodeloutput!(axs, A, t, nquantiles; color, kwargs...)
    return nothing
end

function _plotmodeloutput!(axs, A, t::AbstractVector, nquantiles; kwargs...)
    medquantile = convert(Int, nquantiles / 2 + 0.5)
    for i in medquantile:-1:1
        i >= medquantile && continue 
        _plotmodeloutputband!(
            axs, A, t, i, nquantiles + 1 - i, 2 * i / nquantiles; 
            bandkws(; skip=[:alpha], kwargs...)...
        )
    end
    _plotmodeloutputmedian!(axs, A, t, medquantile; lineskws(; kwargs...)...)
    return nothing
end

function _plotmodeloutputband!(axs, A, t, qlow, qhigh, alpha; kwargs...)
    for (j, ax) in enumerate(axs)
        band!(ax, t, A[:, j, qlow], A[:, j, qhigh]; alpha, kwargs...)
    end
    return nothing
end

function _plotmodeloutputmedian!(axs, A, t, medquantile; kwargs...)
    for (j, ax) in enumerate(axs)
        lines!(ax, t, A[:, j, medquantile]; kwargs...)
    end
    return nothing
end





# Error messages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

