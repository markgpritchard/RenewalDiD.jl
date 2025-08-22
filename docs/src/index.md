
```@meta
DocTestSetup = quote
    using RenewalDiD

    ENV["COLUMNS"] = 100
end
```

# RenewalDiD.jl

Functions to fit parameters using difference-in-differences methods to a renewal equation 
    model. Also includes code to generate simulated datasets.

## Intervention arrays

Easy to generate and view `AbstractArray`s that can be used to indicate intervention times.

```@docs
AbstractInterventionArray
AbstractInterventionVector
InterventionVector
InterventionVector(::AbstractInterventionVector)
AbstractInterventionMatrix
InterventionMatrix
AbstractInterventionArray3
InterventionArray
```

## Generation intervals 

Functions to provide a generation interval for the renewal equation.

```@docs
g_covid
g_seir
generationtime
testgenerationtime
vectorg_seir
```

## Renewal DiD data 

Structures for passing data to the function `renewaldid`.

```@docs
AbstractRenewalDiDData
RenewalDiDData
RenewalDiDDataUnlimitedPopn
```
need to add `SimulationData`

