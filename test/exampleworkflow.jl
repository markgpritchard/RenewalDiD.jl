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
        u0=u0_1, beta=_beta1, gamma, delta, theta, sigma, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=_beta2, gamma, delta, theta, sigma, intervention=30,
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=_beta3, gamma, delta, theta, sigma, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3)
end

model1 = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors(; alphaprior=Normal(log(2.5), 1), sigma_thetaprior=Exponential(0.075));                          
    gamma=0.2, sigma=0.5               
)

priorschain = sample(rng, model1, Prior(), 1000)
priorsdf = DataFrame(priorschain)
priortraceplot = trplot(priorsdf, :tau)
priorsfittedoutputs = samplerenewaldidinfections(
    priorsdf, sim, g_seir; 
    gamma=0.2, sigma=0.5
)
priorsoutputquantiles = quantilerenewaldidinfections(
    priorsfittedoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsplot = plotmodel(priorsoutputquantiles, sim)

initindices = findall(x -> x > 996, ordinalrank(priorsdf.lp)) 
priorsfittedinitoutputs = samplerenewaldidinfections(
    priorsdf, sim, g_seir, initindices; 
    gamma=0.2, sigma=0.5
)
priorsoutputinitquantiles = quantilerenewaldidinfections(
    priorsfittedinitoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
priorsinitplot = plotmodel(priorsoutputinitquantiles, sim)

priorinitparams = [[values(priorsdf[i, 3:733])...] for i in initindices]

shortchain = sample(
    rng, model1, NUTS(0.65; adtype=AutoReverseDiff()), MCMCThreads(), 20, 4; 
    initial_params=priorinitparams,
) 
shortdf = DataFrame(shortchain)
p1 = trplot(shortdf, :tau)
p2 = tracerankplot(shortdf, :tau; binsize=4)

shortfittedoutputs = samplerenewaldidinfections(
    shortdf, sim, g_seir; gamma=0.2, sigma=0.5
)
shortoutputquantiles = quantilerenewaldidinfections(
    shortfittedoutputs, [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
)
p3 = plotmodel(shortoutputquantiles, sim)



#=

chain = sample(rng, model1, NUTS(1000, 0.65), MCMCThreads(), 20, 4;) 
df = DataFrame(chain)

RenewalDiD.Plotting.traceplot(df, :tau)
tracerankplot(df, :tau)

fittedoutputs = samplerenewaldidinfections(g_seir, df; data=sim, gamma=0.2, sigma=0.5)
outputquantiles = quantilerenewaldidinfections(fittedoutputs, [0.025, 0.05, 0.5, 0.95, 0.975])

plotmodel(outputquantiles, sim)
=#

@testset "everything ran and produced expected types" begin
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
