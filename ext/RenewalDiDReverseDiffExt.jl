# vcat with AutoReverseDiff

module RenewalDiDReverseDiffExt 

import RenewalDiD
import ReverseDiff

function RenewalDiD._myvcat(v::AbstractVector{<:Number}, x::ReverseDiff.TrackedReal)
    return ReverseDiff.track(vcat, v..., x)
end

end  # module RenewalDiDReverseDiffExt 
