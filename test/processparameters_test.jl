# test functions for processing output of fitted parameters

using RenewalDiD
using Test

u0_1 = simulationu0(; s=99_990, e=10)
u0_2 = simulationu0(; s=199_990, e=10)
u0_3 = simulationu0(; s=49_993, e=7)
beta1(t) = 0.3 + 0.1 * cos((t - 20) * 2pi / 365)
beta2(t) = beta1(t) * t < 40 ? 1.1 : 0.88
beta3(t) = beta1(t) * t < 60 ? 0.9 : 0.72
s1 = packsimulationtuple( ; 
    u0=u0_1, beta=beta1, gamma=0.2, delta=0.3, theta=0.6, sigma=0.5, intervention=nothing,
)
s2 = packsimulationtuple( ; 
    u0=u0_2, beta=beta2, gamma=0.2, delta=0.3, theta=0.6, sigma=0.5, intervention=40,
)
s3 = packsimulationtuple( ; 
    u0=u0_3, beta=beta3, gamma=0.2, delta=0.3, theta=0.6, sigma=0.5, intervention=60,
)
sim = packsimulations(100, s1, s2, s3)
