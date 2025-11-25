# test that an example workflow throws no errors and returns expected types

import ReverseDiff

using CairoMakie
using Random
using RenewalDiD
using RenewalDiD.Plotting
using StatsBase
using Test
using Turing

rng = Xoshiro(1729)

sim = let 
    u0_1 = simulationu0(; S=100_000, E=5, I=5)
    u0_2 = simulationu0(; S=200_000, E=5, I=5)
    u0_3 = simulationu0(; S=50_000, E=5, I=5)
    eta = 0.2
    phi = 0.6
    sigma = 0.5
    _beta1counter(t) = 0.4 + 0.1 * cos((t - 20) * 2pi / 365) 
    _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t)
    _beta2(t) = t >= 30 ? 0.72 * _beta1counter(t) : 0.9 * _beta1counter(t)
    _beta3(t) = 1.05 * _beta1counter(t)
    s1 = packsimulationtuple( ; 
        u0=u0_1, beta=_beta1, sigma, eta, phi, intervention=[50, 70],
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=_beta2, sigma, eta, phi, intervention=[30, nothing],
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=_beta3, sigma, eta, phi, intervention=[nothing, 40],
    )
    packsimulations(rng, 100, s1, s2, s3; sampletime=14)
end

model1 = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2.5), 1), 
        sigma_gammaprior=Exponential(0.2),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(6, 4),
        delaydistn=LogNormal(log(5), log(2)),
    );                          
    mu=0.2, kappa=0.5               
)

priorschain = sample(rng, model1, Prior(), 1_000; progress=false)
priorsdf = DataFrame(priorschain)
priortraceplot = trplot(priorsdf; ncols=5, nplots=50, size=(1000, 1000))  # examine 50 variables
priorsfittedoutputs = samplerenewaldidinfections(model1, priorsdf; mu=0.2, kappa=0.5,)
priorsoutputquantiles = quantilerenewaldidinfections(
    priorsfittedoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot = plotmodel(
    priorsoutputquantiles, sim; 
    linewidth=1, interventionlinestyle=(:dot, :dense)
)

initindices = findall(x -> x <= 4, ordinalrank(priorsdf.lp; rev=true)) 
priorsfittedinitoutputs = samplerenewaldidinfections(
    model1, priorsdf, initindices; 
    mu=0.2, kappa=0.5,
)
priorsoutputinitquantiles = quantilerenewaldidinfections(
    priorsfittedinitoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsinitplot = plotmodel(priorsoutputinitquantiles, sim)

priorinitparams = [[values(priorsdf[i, 3:430])...] for i in initindices]

map_estimate = maximum_a_posteriori(model1; adtype=AutoReverseDiff(), maxtime=60)
map_df = map_DataFrame(map_estimate)
map_outputs = samplerenewaldidinfections(model1, map_df; mu=0.2, kappa=0.5,)
map_quantiles = quantilerenewaldidinfections(map_outputs, [0.5])
map_plot = plotmodel(map_quantiles, sim)

# do not test the MCMC as this takes a long time and does not add beyond existing tests

@testset "everything ran and produced expected types" begin  # to add further tests later
    @test priorschain isa Chains
    @test priorsdf isa DataFrame 
    @test priortraceplot isa Figure 
    @test priorsplot isa Figure 
    @test map_df isa DataFrame
    @test map_plot isa Figure
end
