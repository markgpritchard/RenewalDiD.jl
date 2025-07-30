# plot output and data from model 

# Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combined figure 

function RenewalDiD.Plotting.plotmodel(args...; kwargs...)
    fig = Figure(; kwargs...)
    plotmodel!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.Plotting.plotmodel!(
    gl::FigOrGridLayout, modeloutputs::AbstractArray, args...; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, modeloutputs; axiskws(; kwargs...)...)
    plotmodel!(axs, modeloutputs, args...; kwargs...)
    linkaxes!(axs...)
    return axs
end

function RenewalDiD.Plotting.plotmodel!(
    axs::AbstractVector{Axis}, 
    modeloutputs::AbstractArray, 
    data::RenewalDiDData, 
    t=automatic; 
    kwargs...
) 
    observedcases = data.observedcases
    interventions = data.interventions
    return plotmodel!(
        axs, modeloutputs, observedcases, interventions, t; 
        kwargs...
    ) 
end

function RenewalDiD.Plotting.plotmodel!(
    axs::AbstractVector{Axis}, 
    modeloutputs::AbstractArray, 
    observedcases::AbstractArray, 
    interventions::AbstractArray, 
    t=automatic; 
    kwargs...
) 
    return _plotmodel!(axs, modeloutputs, observedcases, interventions, t; kwargs...)
end

function _plotmodel!(
    axs, modeloutputs, observedcases, interventions, t; 
    color=automatic,    
    datacolor=automatic,
    interventioncolor=automatic,
    modeloutputcolor=automatic,
    kwargs...
) 
    plotmodeloutput!(
        axs, modeloutputs, t; 
        _plotmodel_plotoutputkws(color, modeloutputcolor; kwargs...)...
    ) 
    plotmodeldata!(
        axs, observedcases, t; 
        _plotmodel_plotdatakws(color, datacolor; kwargs...)...
    ) 
    plotmodelintervention!(
        axs, interventions; 
        _plotmodel_plotinterventionkws(color, interventioncolor; kwargs...)...
    ) 
    return axs
end

_plotmodel_plotoutputkws(color, ::Any; kwargs...) = (color, kwargs...)

function _plotmodel_plotoutputkws(::Automatic, modeloutputcolor; kwargs...)
    return (color=modeloutputcolor, kwargs...)
end

_plotmodel_plotoutputkws(::Automatic, ::Automatic; kwargs...) = (color=Cycled(1), kwargs...)

_plotmodel_plotdatakws(color, ::Any; kwargs...) = scatterkws(; color, kwargs...)

function _plotmodel_plotdatakws(::Automatic, datacolor; kwargs...)
    return scatterkws(; color=datacolor, kwargs...)
end

_plotmodel_plotdatakws(::Automatic, ::Automatic; kwargs...) = scatterkws(; kwargs...)

_plotmodel_plotinterventionkws(color, ::Any; kwargs...) = vlineskws(; color, kwargs...)

function _plotmodel_plotinterventionkws(::Automatic, interventioncolor; kwargs...)
    return vlineskws(; color=interventioncolor, kwargs...)
end

_plotmodel_plotinterventionkws(::Automatic, ::Automatic; kwargs...) = vlineskws(; kwargs...)

_plotmodeloutputaxs(gl, A; kwargs...) = [Axis(gl[1, i]; kwargs...) for i in axes(A, 2)]

## Model output 

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
    return axs
end

function RenewalDiD.Plotting.plotmodeloutput!(
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

function __plotmodeloutput!(axs, A, t, nquantiles; kwargs...)
    medquantile = convert(Int, nquantiles / 2 + 0.5)
    for i in medquantile:-1:1
        i >= medquantile && continue 
        _plotmodeloutputband!(
            axs, A, t, i, nquantiles + 1 - i, 2 * i / nquantiles; 
            bandkws(; skip=[:alpha], kwargs...)...
        )
    end
    _plotmodeloutputmedian!(axs, A, t, medquantile; lineskws(; kwargs...)...)
    return axs
end

function _plotmodeloutputband!(axs, A, t, qlow, qhigh, alpha; kwargs...)
    for (j, ax) in enumerate(axs)
        band!(ax, t, A[:, j, qlow], A[:, j, qhigh]; alpha, kwargs...)
    end
    return axs
end

function _plotmodeloutputmedian!(axs, A, t, medquantile; kwargs...)
    for (j, ax) in enumerate(axs)
        lines!(ax, t, A[:, j, medquantile]; kwargs...)
    end
    return axs
end

## Model data 

function RenewalDiD.Plotting.plotmodeldata(args...; kwargs...)
    fig = Figure(; kwargs...)
    plotmodeldata!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.Plotting.plotmodeldata!(
    gl::FigOrGridLayout, A::AbstractArray, t=automatic; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, A; axiskws(; kwargs...)...)
    _plotmodeldata!(axs, A, t; scatterkws(; kwargs...)...)
    linkaxes!(axs...)
    return axs
end

function RenewalDiD.Plotting.plotmodeldata!(
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

function _plotmodeldatascatter!(axs, A, t; color=:black, marker=:x, markersize=3, kwargs...)
    for (j, ax) in enumerate(axs)
        scatter!(ax, t, A[:, j]; color, marker, markersize, kwargs...)
    end
    return axs
end

## Plot time of intervention 

function RenewalDiD.Plotting.plotmodelintervention(args...; kwargs...)
    fig = Figure(; kwargs...)
    plotmodelintervention!(fig, args...; kwargs...)
    return fig
end

function RenewalDiD.Plotting.plotmodelintervention!(
    gl::FigOrGridLayout, A::AbstractArray; 
    kwargs...
) 
    axs = _plotmodeloutputaxs(gl, A; axiskws(; kwargs...)...)
    plotmodelintervention!(axs, A; vlineskws(; kwargs...)...)
    linkaxes!(axs...)
    return axs
end

function RenewalDiD.Plotting.plotmodelintervention!(
    axs::AbstractVector{Axis}, A::AbstractArray; 
    kwargs...
) 
    length(axs) == size(A, 2) || throw(ArgumentError("to do"))
    return _plotmodelintervention!(axs, A; kwargs...)
end

function _plotmodelintervention!(axs, A; kwargs...)
    return _plotmodelinterventionvlines!(axs, A; kwargs...)
end

function _plotmodelinterventionvlines!(axs, A::InterventionMatrix; color=:red, kwargs...)
    for (j, ax) in enumerate(axs)
        isnothing(RenewalDiD._interventionstarttimes(A, j)) && continue
        vlines!(ax, RenewalDiD._interventionstarttimes(A, j); color, kwargs...)
    end
    return axs
end

