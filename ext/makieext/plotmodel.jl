# plot output and data from model 

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combined figure 

function RenewalDiD.plotmodel(args...; kwargs...)
    fig = Figure(; kwargs...)
    RenewalDiD.plotmodel!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.plotmodel!(
    gl::FigOrGridLayout, modeloutputs::AbstractArray, args...; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, modeloutputs; kwargs...)
    RenewalDiD.plotmodel!(axs, modeloutputs, args...; kwargs...)
    linkaxes!(axs...)
    return axs
end

function RenewalDiD.plotmodel!(
    axs::AbstractVector{Axis}, 
    modeloutputs::AbstractArray, 
    data::RenewalDiDData, 
    t=automatic; 
    kwargs...
) 
    observedcases = data.observedcases
    interventions = data.interventions
    return RenewalDiD.plotmodel!(
        axs, modeloutputs, observedcases, interventions, t; 
        kwargs...
    ) 
end

function RenewalDiD.plotmodel!(
    axs::AbstractVector{Axis}, 
    modeloutputs::AbstractArray, 
    observedcases::AbstractArray, 
    interventions::AbstractArray, 
    t=automatic; 
    kwargs...
) 
    return _plotmodel!(axs, modeloutputs, observedcases, interventions, t; kwargs...)
end

function _plotmodel!(axs, modeloutputs, observedcases, interventions, t; kwargs...) 
    RenewalDiD.plotmodeloutput!(axs, modeloutputs, t; kwargs...)
    RenewalDiD.plotmodeldata!(axs, observedcases, t; kwargs...)
    RenewalDiD.plotmodelintervention!(axs, interventions; kwargs...)
    return axs
end

function _plotmodeloutputaxs(gl, A; kwargs...)
    return [Axis(gl[1, i]; axiskws(; prefix=:axis, kwargs...)...) for i in axes(A, 2)]
end

## Model output 

function RenewalDiD.plotmodeloutput(args...; kwargs...)
    fig = Figure(; kwargs...)
    RenewalDiD.plotmodeloutput!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.plotmodeloutput!(
    gl::FigOrGridLayout, A::AbstractArray, t=automatic; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, A; kwargs...)
    _plotmodeloutput!(axs, A, t; kwargs...)
    linkaxes!(axs...)
    return axs
end

function RenewalDiD.plotmodeloutput!(
    axs::AbstractVector{Axis}, A::AbstractArray, t=automatic; 
    kwargs...
) 
    length(axs) == size(A, 2) || throw(ArgumentError("to do"))
    return _plotmodeloutput!(axs, A, t; kwargs...)
end

function _plotmodeloutput!(axs, A, ::Automatic; kwargs...)
    t = axes(A, 1)
    return _plotmodeloutput!(axs, A, t; kwargs...)
end

function _plotmodeloutput!(axs, A, t::AbstractVector; color=Cycled(1), kwargs...)
    # set colour so the median line and credible interval bands have the same colour
    return __plotmodeloutput!(axs, A, t; color, kwargs...)
end

function __plotmodeloutput!(axs, A::AbstractMatrix, t; kwargs...) 
    nquantiles = 1  # Matrix so only one dataset to plot (i.e. median only)
    return __plotmodeloutput!(axs, A, t, nquantiles; kwargs...)
end

function __plotmodeloutput!(axs, A::AbstractArray{T, 3}, t; kwargs...) where T
    nquantiles = size(A, 3)
    isodd(size(A, 3)) || throw(ArgumentError("to do"))
    return __plotmodeloutput!(axs, A, t, nquantiles; kwargs...)
end

function __plotmodeloutput!(axs, A, t, nquantiles; modelcolor=Cycled(1), kwargs...)
    medquantile = convert(Int, nquantiles / 2 + 0.5)
    for i in medquantile:-1:1
        i >= medquantile && continue 
        _plotmodeloutputband!(
            axs, 
            A, 
            t, 
            i,  # qlow
            nquantiles + 1 - i,  # qhigh
            2 * i / nquantiles;  # alpha calculated here, not passed as keyword argument
            modelcolor, kwargs...
        )
    end
    _plotmodeloutputmedian!(axs, A, t, medquantile; modelcolor, kwargs...)
    return axs
end

function _plotmodeloutputband!(axs, A, t, qlow, qhigh, alpha; kwargs...)
    for (j, ax) in enumerate(axs)
        band!(
            ax, t, A[:, j, qlow], A[:, j, qhigh]; 
            alpha, bandkws(; prefix=:model, skip=[:alpha], kwargs...)...
        )
    end
    return axs
end

function _plotmodeloutputmedian!(axs, A, t, medquantile; kwargs...)
    for (j, ax) in enumerate(axs)
        lines!(ax, t, A[:, j, medquantile]; lineskws(; prefix=:model, kwargs...)...)
    end
    return axs
end

## Model data 

function RenewalDiD.plotmodeldata(args...; kwargs...)
    fig = Figure(; kwargs...)
    RenewalDiD.plotmodeldata!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.plotmodeldata!(
    gl::FigOrGridLayout, A::AbstractArray, t=automatic; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, A; kwargs...)
    _plotmodeldata!(axs, A, t; kwargs...)
    linkaxes!(axs...)
    return axs
end

function RenewalDiD.plotmodeldata!(
    axs::AbstractVector{Axis}, A::AbstractArray, t=automatic; 
    kwargs...
) 
    length(axs) == size(A, 2) || throw(ArgumentError("to do"))
    return _plotmodeldata!(axs, A, t; kwargs...)
end

function _plotmodeldata!(axs, A, ::Automatic; kwargs...)
    t = axes(A, 1)
    return _plotmodeldata!(axs, A, t; kwargs...)
end

function _plotmodeldata!(axs, A, t::AbstractVector; kwargs...)
    return _plotmodeldatascatter!(axs, A, t; kwargs...)
end

function _plotmodeldatascatter!(
    axs, A, t; 
    datacolor=:black, datamarker=:x, datamarkersize=3, kwargs...
)
    for (j, ax) in enumerate(axs)
        scatter!(
            ax, t, A[:, j]; 
            scatterkws(; prefix=:data, datacolor, datamarker, datamarkersize, kwargs...)...
        )
    end
    return axs
end

## Plot time of intervention 

function RenewalDiD.plotmodelintervention(args...; kwargs...)
    fig = Figure(; kwargs...)
    RenewalDiD.plotmodelintervention!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.plotmodelintervention!(
    gl::FigOrGridLayout, A::AbstractArray; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, A; kwargs...)
    RenewalDiD.plotmodelintervention!(axs, A; kwargs...)
    linkaxes!(axs...)
    return axs
end

function RenewalDiD.plotmodelintervention!(
    axs::AbstractVector{Axis}, A::AbstractArray; 
    kwargs...
) 
    length(axs) == size(A, 2) || throw(ArgumentError("to do"))
    return _plotmodelintervention!(axs, A; kwargs...)
end

function _plotmodelintervention!(axs, A; kwargs...)
    return _plotmodelinterventionvlines!(axs, A; kwargs...)
end

function _plotmodelinterventionvlines!(
    axs, A::AbstractMatrix; 
    interventioncolor=:red, kwargs...
)
    for (j, ax) in enumerate(axs)
        isnothing(RenewalDiD._interventionstarttimes(A, j)) && continue
        vlines!(
            ax, RenewalDiD._interventionstarttimes(A, j); 
            vlineskws(; prefix=:intervention, interventioncolor, kwargs...)...
        )
    end
    return axs
end

