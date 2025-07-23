# The package includes functions for plotting outputs in CairoMakie. To avoid loading this
# package unnecessarily for users who prefer to use other plotting packages, this file
# simply contains the name of each function to be exported. Methods for these functions are
# exported only if the user is also `using CairoMakie`

function traceplot end
function traceplot! end
function tracerankplot end
function tracerankplot! end
