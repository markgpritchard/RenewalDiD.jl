# test the function `renewaldid`

import ReverseDiff  # `Mooncake` performs very slowly in this set of tests

using DynamicPPL.TestUtils.AD: run_ad
using Random: Xoshiro
using RenewalDiD
using Test
using Turing: AutoReverseDiff, Beta, Exponential, LogNormal, Normal, Prior, sample

rng = Xoshiro(1729)

sim = let 
    u0_1 = simulationu0(; S=100_000, E=5, I=5)
    u0_2 = simulationu0(; S=200_000, E=5, I=5)
    u0_3 = simulationu0(; S=50_000, E=5, I=5)
    eta = 0.2
    phi = 0.6
    sigma = 0.5
    _beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365) 
    _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t)
    _beta2(t) = t >= 30 ? 0.72 * _beta1counter(t) : 0.9 * _beta1counter(t)
    _beta3(t) = 1.05 * _beta1counter(t)
    s1 = packsimulationtuple(; u0=u0_1, beta=_beta1, eta, phi, sigma, intervention=50,)
    s2 = packsimulationtuple(; u0=u0_2, beta=_beta2, eta, phi, sigma, intervention=30,)
    s3 = packsimulationtuple(; u0=u0_3, beta=_beta3, eta, phi, sigma, intervention=nothing,)
    packsimulations(rng, 100, s1, s2, s3; sampletime=14)
end

sim2 = let  # unlimited population
    u0_1 = simulationu0(; S=100_000, E=5, I=5)
    u0_2 = simulationu0(; S=200_000, E=5, I=5)
    u0_3 = simulationu0(; S=50_000, E=5, I=5)
    eta = 0.2
    phi = 0.6
    sigma = 0.5
    _beta1counter(t) = 0.4 + 0.1 * cos((t-20) * 2pi / 365) 
    _beta1(t) = t >= 50 ? 0.8 * _beta1counter(t) : _beta1counter(t)
    _beta2(t) = t >= 30 ? 0.72 * _beta1counter(t) : 0.9 * _beta1counter(t)
    _beta3(t) = 1.05 * _beta1counter(t)
    s1 = packsimulationtuple(; u0=u0_1, beta=_beta1, eta, phi, sigma, intervention=50,)
    s2 = packsimulationtuple(; u0=u0_2, beta=_beta2, eta, phi, sigma, intervention=30,)
    s3 = packsimulationtuple(; u0=u0_3, beta=_beta3, eta, phi, sigma, intervention=nothing,)
    packsimulations(rng, 100, s1, s2, s3; sampletime=14, unlimitedpop=true)
end

model1 = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2.5), 1), 
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(6, 4),
        delaydistn=LogNormal(log(2.5), log(5)),
    );                          
    mu=0.2, kappa=0.5,               
)

model2 = renewaldid(                      
    sim2, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2.5), 1), 
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(6, 4),
        delaydistn=LogNormal(log(2.5), log(5)),
    );                          
    mu=0.2, kappa=0.5,               
)

model3 = renewaldid(                      
    sim, 
    g_seir, 
    RenewalDiDPriors( ; 
        alphaprior=Normal(log(2.5), 1), 
        sigma_thetaprior=Exponential(0.075), 
        psiprior=Beta(6, 4),
        delaydistn=Exponential(1 / 0.3),
    );                          
    mu=0.2, kappa=0.5,               
)

@testset "sufficient non-`NaN` log likelihoods from prior" begin
    # test would have been failed by previous version of function `_renewaldid`
    ps = [sample(rng, m, Prior(), 1000; progress=false) for m in [model1, model2, model3]]
    ds = DataFrame.(ps)
    for d in ds 
        nonnanloglikelihood = [isnan(x) ? 0 : 1 for x in d.loglikelihood]
        @test sum(nonnanloglikelihood) > 500 
    end
end

@testset "any `NaN` gradients in model 1?" begin
    # this test seems dependent on the rng supplied -- would be good to clarify why and make
    # more robust
    adtype = AutoReverseDiff()
    result = run_ad(model1, adtype; rng=Xoshiro(2000), test=false, verbose=false,);
    @test sum(isnan.(result.grad_actual)) == 0
    @test isnothing(findfirst(isnan, result.grad_actual))
end

@testset "any `NaN` gradients in model 2?" begin
    adtype = AutoReverseDiff()
    result = run_ad(model2, adtype; test=false, verbose=false,);
    @test sum(isnan.(result.grad_actual)) == 0
    @test isnothing(findfirst(isnan, result.grad_actual))
end

@testset "any `NaN` gradients in model 3?" begin
    adtype = AutoReverseDiff()
    result = run_ad(model3, adtype; rng=Xoshiro(2000), test=false, verbose=false,);
    @test sum(isnan.(result.grad_actual)) == 0
    @test isnothing(findfirst(isnan, result.grad_actual))
end

@testset "inferred type" begin
    @test @inferred model1.f(model1.args...) == model1.f(model1.args...)
end
