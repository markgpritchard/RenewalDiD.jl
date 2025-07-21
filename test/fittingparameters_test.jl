# test the function `renewaldid`

using RenewalDiD
using Test
using StableRNGs
using Turing

# initially just a test that it doesn't error

rng = StableRNG(1)

sim = let 
    u0_1 = simulationu0(; s=98, e=2)
    u0_2 = simulationu0(; s=198, e=2)
    u0_3 = simulationu0(; s=48, e=2)
    beta1(t) = 0.3 + 0.1 * cos((t - 20) * 2pi / 365)
    beta2(t) = beta1(t) * t < 4 ? 1.1 : 0.88
    beta3(t) = beta1(t) * t < 6 ? 0.9 : 0.72
    gamma = 0.2
    delta = 0.3
    theta = 0.6
    sigma = 0.5
    s1 = packsimulationtuple( ; 
        u0=u0_1, beta=beta1, gamma, delta, theta, sigma, intervention=nothing,
    )
    s2 = packsimulationtuple( ; 
        u0=u0_2, beta=beta2, gamma, delta, theta, sigma, intervention=4,
    )
    s3 = packsimulationtuple( ; 
        u0=u0_3, beta=beta3, gamma, delta, theta, sigma, intervention=6,
    )
    packsimulations(rng, 10, s1, s2, s3)
end

model =  renewaldid(
    sim, g_seir, packpriors(; sigma_thetaprior=Exponential(0.05)); 
    gamma=0.2, sigma=0.5
)

chain = sample(rng, model, NUTS(), MCMCThreads(), 20, 4; verbose=false, progress=false)

df = DataFrame(chain)

@test df.tau[1] == -0.04157675192947863
