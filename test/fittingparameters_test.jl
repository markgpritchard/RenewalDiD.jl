# test the function `renewaldid`

using DynamicPPL.TestUtils.AD: run_ad
using Random
using RenewalDiD
using ReverseDiff
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
    _beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365) 
    _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t)
    _beta2(t) = t >= 30 ? 0.72 * _beta1counter(t) : 0.9 * _beta1counter(t)
    _beta3(t) = 1.05 * _beta1counter(t)
    s1 = packsimulationtuple( ; 
        u0=u0_1, beta=_beta1, mu, delta, psi, kappa, intervention=50,
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=_beta2, mu, delta, psi, kappa, intervention=30,
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=_beta3, mu, delta, psi, kappa, intervention=nothing,
    )
    packsimulations(rng, 100, s1, s2, s3; sampletime=14)
end

model1 = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2.5), 1), 
        mu_delayprior=log(2.5),  # cannot currently accept distributions
        sigma_delayprior=log(5),  # cannot currently accept distributions
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(6, 4)
    );                          
    mu=0.2, kappa=0.5               
)

@testset "any `NaN` gradients?" begin
    adtype = AutoReverseDiff()
    result = run_ad(model1, adtype; test=false, verbose=false,);
    @test sum(isnan.(result.grad_actual)) == 0
    @test isnothing(findfirst(isnan, result.grad_actual))
end

@testset "mode estimate" begin
    _t = time()
    map_estimate = maximum_likelihood(model1; adtype=AutoReverseDiff(), maxtime=30)
    @test time() - _t < 50 
    map_df = map_DataFrame(map_estimate)
    @test map_df isa DataFrame
end
