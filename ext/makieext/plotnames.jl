# outputs from MCMC that are not plotted 

const NOPLOTNAMES = [ 
    "iteration", 
    "chain", 
    "lp", 
    "n_steps", 
    "is_accept", 
    "acceptance_rate", 
    "log_density", 
    "hamiltonian_energy", 
    "hamiltonian_energy_error", 
    "max_hamiltonian_energy_error", 
    "tree_depth", 
    "numerical_error", 
    "step_size", 
    "nom_step_size"
]

function _plotnames(df) 
    allnames = names(df)
    plotnames = allnames[findall(x -> x ∉ NOPLOTNAMES, allnames)]
    return plotnames 
end

_plotnames(df, vs) = _plotnames(df)[vs]
