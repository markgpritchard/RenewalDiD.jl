# assign appropriate keyword arguments to plotting functions 

# Functions in Makie.jl accept similar keyword arguments but will throw an error if 
# inappropriate arguments are provided.  The functions defined in this extension pass 
# keyword arguments to multiple functions.  Functions on this page provide the correct 
# keyword argument for each function.

# Nothing in this file is exported

# is this a good way of doing it?

const AXISKWS = [
    :alignmode
    :aspect
    :autolimitaspect    
    :backgroundcolor    
    :bottomspinecolor   
    :bottomspinevisible 
    :dim1_conversion    
    :dim2_conversion    
    :flip_ylabel        
    :halign
    :height
    :leftspinecolor     
    :leftspinevisible   
    :limits
    :palette
    :panbutton
    :rightspinecolor    
    :rightspinevisible  
    :spinewidth
    :subtitle
    :subtitlecolor      
    :subtitlefont       
    :subtitlegap        
    :subtitlelineheight 
    :subtitlesize       
    :subtitlevisible    
    :tellheight
    :tellwidth
    :title
    :titlealign
    :titlecolor
    :titlefont
    :titlegap
    :titlelineheight   
    :titlesize
    :titlevisible       
    :topspinecolor      
    :topspinevisible    
    :valign
    :width
    :xautolimitmargin   
    :xaxisposition      
    :xgridcolor
    :xgridstyle
    :xgridvisible       
    :xgridwidth
    :xlabel
    :xlabelcolor        
    :xlabelfont
    :xlabelpadding      
    :xlabelrotation     
    :xlabelsize
    :xlabelvisible      
    :xminorgridcolor    
    :xminorgridstyle    
    :xminorgridvisible  
    :xminorgridwidth    
    :xminortickalign    
    :xminortickcolor    
    :xminorticks        
    :xminorticksize     
    :xminorticksvisible 
    :xminortickwidth    
    :xpankey
    :xpanlock
    :xrectzoom
    :xreversed
    :xscale
    :xtickalign
    :xtickcolor
    :xtickformat        
    :xticklabelalign    
    :xticklabelcolor    
    :xticklabelfont     
    :xticklabelpad      
    :xticklabelrotation 
    :xticklabelsize     
    :xticklabelspace    
    :xticklabelsvisible 
    :xticks
    :xticksize
    :xticksmirrored     
    :xticksvisible      
    :xtickwidth
    :xtrimspine
    :xzoomkey
    :xzoomlock
    :yautolimitmargin   
    :yaxisposition      
    :ygridcolor
    :ygridstyle
    :ygridvisible       
    :ygridwidth
    :ylabel
    :ylabelcolor        
    :ylabelfont
    :ylabelpadding      
    :ylabelrotation     
    :ylabelsize
    :ylabelvisible      
    :yminorgridcolor    
    :yminorgridstyle    
    :yminorgridvisible  
    :yminorgridwidth    
    :yminortickalign    
    :yminortickcolor    
    :yminorticks        
    :yminorticksize     
    :yminorticksvisible 
    :yminortickwidth    
    :ypankey
    :ypanlock
    :yrectzoom
    :yreversed
    :yscale
    :ytickalign
    :ytickcolor
    :ytickformat        
    :yticklabelalign    
    :yticklabelcolor    
    :yticklabelfont     
    :yticklabelpad      
    :yticklabelrotation 
    :yticklabelsize     
    :yticklabelspace    
    :yticklabelsvisible 
    :yticks
    :yticksize
    :yticksmirrored     
    :yticksvisible      
    :ytickwidth
    :ytrimspine
    :yzoomkey
    :yzoomlock
    :zoombutton
]

const BANDKWS = [
    :alpha
    :backlight
    :clip_planes
    :color
    :colormap
    :colorrange
    :colorscale
    :cycle
    :depth_shift
    :diffuse
    :direction
    :fxaa
    :highclip
    :inspectable
    :inspector_clear
    :inspector_hover
    :inspector_label
    :interpolate
    :lowclip
    :matcap
    :material
    :model
    :nan_color
    :overdraw
    :shading
    :shininess
    :space
    :specular
    :ssao
    :transformation
    :transparency
    :uv_transform
    :visible
]

const LINEKWS = [
    :alpha
    :clip_planes
    :color
    :colormap
    :colorrange
    :colorscale
    :cycle
    :depth_shift
    :fxaa
    :highclip
    :inspectable
    :inspector_clear
    :inspector_hover
    :inspector_label
    :joinstyle
    :linecap
    :linestyle
    :linewidth
    :lowclip
    :miter_limit
    :model
    :nan_color
    :overdraw
    :space
    :ssao
    :transformation
    :transparency
    :visible
]

const SCATTERKWS = [
    :alpha
    :clip_planes        
    :color
    :colormap
    :colorrange
    :colorscale
    :cycle
    :depth_shift        
    :depthsorting       
    :distancefield      
    :font
    :fxaa
    :glowcolor
    :glowwidth
    :highclip
    :inspectable        
    :inspector_clear    
    :inspector_hover    
    :inspector_label    
    :lowclip
    :marker
    :marker_offset      
    :markersize
    :markerspace        
    :model
    :nan_color
    :overdraw
    :rotation
    :space
    :ssao
    :strokecolor        
    :strokewidth        
    :transform_marker   
    :transformation     
    :transparency       
    :uv_offset_width    
    :visible
]

const VLINESKWS = [
    :alpha
    :clip_planes
    :color
    :colormap
    :colorrange
    :colorscale
    :cycle
    :depth_shift
    :fxaa
    :highclip
    :inspectable
    :inspector_clear     
    :inspector_hover     
    :inspector_label     
    :linecap
    :linestyle
    :linewidth
    :lowclip
    :model
    :nan_color
    :overdraw
    :space
    :ssao
    :transformation      
    :transparency        
    :visible
    :ymax
    :ymin
]

axiskws(; kwargs...) = _selectkws(AXISKWS; kwargs...)
bandkws(; kwargs...) = _selectkws(BANDKWS; kwargs...)
lineskws(; kwargs...) = _selectkws(LINEKWS; kwargs...)
scatterkws(; kwargs...) = _selectkws(SCATTERKWS; kwargs...)
vlineskws(; kwargs...) = _selectkws(VLINESKWS; kwargs...)

function _selectkws(expectedarguments; skip=Symbol[], kwargs...)
    kwkeys = Symbol[]
    kwvals = Any[]
    for (k, v) in kwargs
        if k in expectedarguments && k ∉ skip
            push!(kwkeys, k)
            push!(kwvals, v)
        end
    end
    return (; (kwkeys .=> kwvals)...)
end
