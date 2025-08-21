
```@meta
DocTestSetup = quote
    using RenewalDiD

    ENV["COLUMNS"] = 100
end
```

# RenewalDiD.jl

Functions to fit parameters using difference-in-differences methods to a renewal equation 
    model. Also includes code to generate simulated datasets.

```@docs
renewaldid
```

## Intervention arrays

```@docs
AbstractInterventionArray{T, N} <: AbstractArray{T, N}

```