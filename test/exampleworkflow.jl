# test that an example workflow throws no errors and returns expected types

using CairoMakie
using Random
using RenewalDiD
using RenewalDiD.Plotting
using ReverseDiff
using StatsBase
using Test
using Turing

rng = Xoshiro(1729)

sim = let 
    u0_1 = simulationu0(; s=100_000, e=5, i_n=3, i_f=2)
    u0_2 = simulationu0(; s=200_000, e=5, i_n=3, i_f=2)
    u0_3 = simulationu0(; s=50_000, e=5, i_n=3, i_f=2)
    mu = 0.2
    delta = 0.3
    psi = 0.6
    kappa = 0.5
    _beta1counter(t) = 0.4 + 0.1 * cos((t - 20) * 2pi / 365) 
    _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t)
    _beta2(t) = t >= 30 ? 0.72 * _beta1counter(t) : 0.9 * _beta1counter(t)
    _beta3(t) = 1.05 * _beta1counter(t)
    s1 = packsimulationtuple( ; 
        u0=u0_1, beta=_beta1, mu, delta, psi, kappa, intervention=[50, 70],
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=_beta2, mu, delta, psi, kappa, intervention=[30, nothing],
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=_beta3, mu, delta, psi, kappa, intervention=[nothing, 40],
    )
    packsimulations(rng, 100, s1, s2, s3; sampletime=14)
end

model1 = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2.5), 1), 
        mu_delayprior=log(5),
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(6, 4)
    );                          
    mu=0.2, kappa=0.5               
)

priorschain = sample(rng, model1, Prior(), 1_000)
priorsdf = DataFrame(priorschain)
priortraceplot = trplot(priorsdf; ncols=5, nplots=50, size=(1000, 1000))  # examine 50 variables
priorsfittedoutputs = samplerenewaldidinfections(
    g_seir, priorsdf, sim; 
    mu=0.2, kappa=0.5,
)
priorsoutputquantiles = quantilerenewaldidinfections(
    priorsfittedoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot = plotmodel(
    priorsoutputquantiles, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)

initindices = findall(x -> x <= 4, ordinalrank(priorsdf.lp; rev=true)) 
priorsfittedinitoutputs = samplerenewaldidinfections(
    g_seir, priorsdf, sim, initindices; 
    mu=0.2, kappa=0.5,
)
priorsoutputinitquantiles = quantilerenewaldidinfections(
    priorsfittedinitoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsinitplot = plotmodel(priorsoutputinitquantiles, sim)

priorinitparams = [[values(priorsdf[i, 3:430])...] for i in initindices]

map_estimate = maximum_likelihood(model1; adtype=AutoReverseDiff(), maxiters=1000,)
map_df = map_DataFrame(map_estimate)
map_outputs = samplerenewaldidinfections(
    g_seir, map_df, sim; 
    mu=0.2, kappa=0.5,
)
map_quantiles = quantilerenewaldidinfections(map_outputs, [0.5])
map_plot = plotmodel(map_quantiles, sim)


shortchain = sample(
    rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 25, 4; 
    #rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 20, 4; 
    #rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 5, 4; 
    #initial_params=priorinitparams,
    initial_params=repeat([map_estimate.values.array]; outer=4)
) 
shortdf = DataFrame(shortchain)
p1 = trplot(shortdf; ncols=5, nplots=50, size=(1000, 1000))  # examine 50 variables
p2 = tracerankplot(shortdf; binsize=5, ncols=5, nplots=50, size=(1000, 1000)) 

shortfittedoutputs = samplerenewaldidinfections(
    g_seir, shortdf, sim; 
    mu=0.2, kappa=0.5,
)
shortoutputquantiles = quantilerenewaldidinfections(
    shortfittedoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
p3 = plotmodel(shortoutputquantiles, sim)

@testset "everything ran and produced expected types" begin  # to add further tests later
    @test rng isa Xoshiro
    @test sim isa RenewalDiDData{Int64, InterventionArray{Int64}} 
    @test shortchain isa Chains
    @test shortdf isa DataFrame 
    @test p1 isa Figure
    @test p2 isa Figure
    @test shortfittedoutputs isa Array{Float64, 3}
    @test shortoutputquantiles isa Array{Float64, 3} 
    @test p3 isa Figure
end
