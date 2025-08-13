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
    gamma = 0.2
    delta = 0.3
    theta = 0.6
    sigma = 0.5
    _beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365) 
    _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t)
    _beta2(t) = t >= 30 ? 0.72 * _beta1counter(t) : 0.9 * _beta1counter(t)
    _beta3(t) = 1.05 * _beta1counter(t)
    s1 = packsimulationtuple( ; 
        u0=u0_1, beta=_beta1, gamma, delta, theta, sigma, intervention=[50, 70],
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=_beta2, gamma, delta, theta, sigma, intervention=[30, 50],
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=_beta3, gamma, delta, theta, sigma, intervention=[nothing, nothing],
    )
    packsimulations(rng, 100, s1, s2, s3; sampletime=14)
end

model1 = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2.5), 1), 
        mu_delayprior=log(5),
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(6, 4)
    );                          
    gamma=0.2, sigma=0.5               
)

priorschain = sample(rng, model1, Prior(), 10_000)
priorsdf = DataFrame(priorschain)
priortraceplot = trplot(priorsdf; ncols=5, nplots=50, size=(1000, 1000))  # examine 50 variables
priorsfittedoutputs = samplerenewaldidinfections(
    g_seir, priorsdf, sim; 
    gamma=0.2, sigma=0.5,
)
priorsoutputquantiles = quantilerenewaldidinfections(
    priorsfittedoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot = plotmodel(priorsoutputquantiles, sim; linewidth=1)

initindices = findall(x -> x <= 4, ordinalrank(priorsdf.lp; rev=true)) 
priorsfittedinitoutputs = samplerenewaldidinfections(
    g_seir, priorsdf, sim, initindices; 
    gamma=0.2, sigma=0.5,
)
priorsoutputinitquantiles = quantilerenewaldidinfections(
    priorsfittedinitoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsinitplot = plotmodel(priorsoutputinitquantiles, sim)

#priorinitparams = [[values(priorsdf[i, 3:735])...] for i in initindices]
priorinitparams = [[values(priorsdf[i, 3:733])...] for i in initindices]

shortchain = sample(
    rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 25, 4; 
    #rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 20, 4; 
    #rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 5, 4; 
    initial_params=priorinitparams,
) 
shortdf = DataFrame(shortchain)
p1 = trplot(shortdf; ncols=5, nplots=50, size=(1000, 1000))  # examine 50 variables
p2 = tracerankplot(shortdf; binsize=5, ncols=5, nplots=50, size=(1000, 1000)) 

shortfittedoutputs = samplerenewaldidinfections(
    g_seir, shortdf, sim; 
    gamma=0.2, sigma=0.5,
)
shortoutputquantiles = quantilerenewaldidinfections(
    shortfittedoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
p3 = plotmodel(shortoutputquantiles, sim)

@testset "everything ran and produced expected types" begin  # to add further tests later
    @test rng isa Xoshiro
    @test sim isa RenewalDiDData{Int64, InterventionMatrix{Int64}} 
    @test shortchain isa Chains
    @test shortdf isa DataFrame 
    @test p1 isa Figure
    @test p2 isa Figure
    @test shortfittedoutputs isa Array{Float64, 3}
    @test shortoutputquantiles isa Array{Float64, 3} 
    @test p3 isa Figure
end
