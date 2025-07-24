# functions for plotting 

module Plotting 

using RenewalDiD

# to avoid loading the `CairoMakie` package for users who are not planning to use it, this 
# module simply exports the names of the plotting functions. Their signatures are added in
# `RenewalDiDCairoMakieExt` after `using CairoMakie`

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            plottinghint = "\nHINT: Plotting functions from `RenewalDiD.Plotting` are only \
                available after `using CairoMakie`"
            plottingfunctionslist = [
                traceplot, traceplot!, 
                tracerankplot, tracerankplot!
            ]

            if exc.f in plottingfunctionslist && 
                isnothing(Base.get_extension(RenewalDiD, :RenewalDiDCairoMakieExt))

                print(io, plottinghint)
            end
        end
    end
    return nothing
end

export traceplot, traceplot!
export tracerankplot, tracerankplot!

function traceplot end
function traceplot! end
function tracerankplot end
function tracerankplot! end

end  # module Plotting
