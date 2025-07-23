# test functions used in other tests

using RenewalDiD.FittedParameterTestFunctions
using StableRNGs
using Test

rng1 = StableRNG(1)
rng2 = StableRNG(1)
rng3 = StableRNG(1)
rng4 = StableRNG(1)

olddf1 = let
    nrows = 2
    DataFrame(
        :iteration => 1:nrows, 
        :chain => ones(Int, nrows), 
        :tau => [0, rand(rng1)], 
        :alpha => [0, rand(rng1)], 
        :sigma_gamma => [0, rand(rng1)], 
        Symbol("gammas_raw[1]") => rand(rng1, nrows), 
        Symbol("gammas_raw[2]") => rand(rng1, nrows), 
        Symbol("thetas_raw[1]") => rand(rng1, nrows), 
        Symbol("thetas_raw[2]") => rand(rng1, nrows),  
        Symbol("thetas_raw[3]") => rand(rng1, nrows), 
        Symbol("thetas_raw[4]") => rand(rng1, nrows), 
        Symbol("thetas_raw[5]") => rand(rng1, nrows), 
        Symbol("thetas_raw[6]") => rand(rng1, nrows), 
        Symbol("thetas_raw[7]") => rand(rng1, nrows), 
        Symbol("thetas_raw[8]") => rand(rng1, nrows),  
        Symbol("thetas_raw[9]") => rand(rng1, nrows), 
        Symbol("sigma_theta") => [0, rand(rng1)], 
        Symbol("M_x[1, 1]") => rand(rng1, nrows), 
        Symbol("M_x[2, 1]") => rand(rng1, nrows), 
        Symbol("M_x[3, 1]") => rand(rng1, nrows),  
        Symbol("M_x[4, 1]") => rand(rng1, nrows),  
        Symbol("M_x[5, 1]") => rand(rng1, nrows), 
        Symbol("M_x[6, 1]") => rand(rng1, nrows),  
        Symbol("M_x[7, 1]") => rand(rng1, nrows), 
        Symbol("M_x[8, 1]") => rand(rng1, nrows), 
        Symbol("M_x[9, 1]") => rand(rng1, nrows),  
        Symbol("M_x[10, 1]") => rand(rng1, nrows), 
        Symbol("M_x[11, 1]") => rand(rng1, nrows), 
        Symbol("M_x[12, 1]") => rand(rng1, nrows), 
        Symbol("M_x[13, 1]") => rand(rng1, nrows),  
        Symbol("M_x[14, 1]") => rand(rng1, nrows),  
        Symbol("M_x[15, 1]") => rand(rng1, nrows),  
        Symbol("M_x[16, 1]") => rand(rng1, nrows),  
        Symbol("M_x[17, 1]") => rand(rng1, nrows), 
        Symbol("M_x[1, 2]") => rand(rng1, nrows), 
        Symbol("M_x[2, 2]") => rand(rng1, nrows), 
        Symbol("M_x[3, 2]") => rand(rng1, nrows), 
        Symbol("M_x[4, 2]") => rand(rng1, nrows),  
        Symbol("M_x[5, 2]") => rand(rng1, nrows), 
        Symbol("M_x[6, 2]") => rand(rng1, nrows), 
        Symbol("M_x[7, 2]") => rand(rng1, nrows),  
        Symbol("M_x[8, 2]") => rand(rng1, nrows),  
        Symbol("M_x[9, 2]") => rand(rng1, nrows),  
        Symbol("M_x[10, 2]") => rand(rng1, nrows),  
        Symbol("M_x[11, 2]") => rand(rng1, nrows),  
        Symbol("M_x[12, 2]") => rand(rng1, nrows), 
        Symbol("M_x[13, 2]") => rand(rng1, nrows),  
        Symbol("M_x[14, 2]") => rand(rng1, nrows), 
        Symbol("M_x[15, 2]") => rand(rng1, nrows),  
        Symbol("M_x[16, 2]") => rand(rng1, nrows),  
        Symbol("M_x[17, 2]") => rand(rng1, nrows), 
        Symbol("M_x[1, 3]") => rand(rng1, nrows), 
        Symbol("M_x[2, 3]") => rand(rng1, nrows), 
        Symbol("M_x[3, 3]") => rand(rng1, nrows), 
        Symbol("M_x[4, 3]") => rand(rng1, nrows), 
        Symbol("M_x[5, 3]") => rand(rng1, nrows), 
        Symbol("M_x[6, 3]") => rand(rng1, nrows),  
        Symbol("M_x[7, 3]") => rand(rng1, nrows), 
        Symbol("M_x[8, 3]") => rand(rng1, nrows),  
        Symbol("M_x[9, 3]") => rand(rng1, nrows),  
        Symbol("M_x[10, 3]") => rand(rng1, nrows),  
        Symbol("M_x[11, 3]") => rand(rng1, nrows), 
        Symbol("M_x[12, 3]") => rand(rng1, nrows),  
        Symbol("M_x[13, 3]") => rand(rng1, nrows),  
        Symbol("M_x[14, 3]") => rand(rng1, nrows),  
        Symbol("M_x[15, 3]") => rand(rng1, nrows),  
        Symbol("M_x[16, 3]") => rand(rng1, nrows), 
        Symbol("M_x[17, 3]") => rand(rng1, nrows),  
        :lp => -100 .* rand(rng1, nrows), 
        :n_steps => ones(nrows),  
        :is_accept => ones(nrows), 
        :acceptance_rate => rand(rng1, nrows), 
        :log_density => -100 .* rand(rng1, nrows),  
        :hamiltonian_energy => rand(rng1, nrows), 
        :hamiltonian_energy_error => rand(rng1, nrows), 
        :max_hamiltonian_energy_error => rand(rng1, nrows), 
        :tree_depth => 3 .* ones(nrows), 
        :numerical_error => zeros(nrows), 
        :step_size => rand(rng1, nrows), 
        :nom_step_size => rand(rng1, nrows),
    )
end

sim1 = let 
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
    packsimulations(rng3, 10, s1, s2, s3)
end

df1 = testdataframe(
    rng2; 
    nchains=1, 
    niterations=2, 
    ngroups=3, 
    ntimes=10, 
    nseeds=7, 
    tau=[0, rand(rng2)], 
    alpha=[0, rand(rng2)], 
    sigma_gamma=[0, rand(rng2)], 
    sigma_theta=[0, rand(rng2)], 
)

@testset for (i, name) in enumerate(names(df1))
    @test name == names(olddf1)[i]
    @testset for j in 1:2 
        if getproperty(df1, name)[j] == round(Int, getproperty(df1, name)[j])
            @test getproperty(df1, name)[j] == getproperty(olddf1, name)[j]
        end
    end
end

@testset "Test simulation" begin
    @test testsimulation(rng4) == sim1
end
