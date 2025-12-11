# test functions that generate simulations 

import Random

using RenewalDiD
using StableRNGs: StableRNG
using Test

function testseirrates(u, t; beta=0, sigma=0, eta=0, phi=0)
    return RenewalDiD._seirrates(u, t, beta, sigma, eta, phi)
end

function testsimulateday!(rng, u, t; beta=0, sigma=0, eta=0, phi=0)
    return RenewalDiD._simulateday!(rng, u, t, beta, sigma, eta, phi)
end

expectedv1 = 0.0028200566597855305
S1 = runsimulation(StableRNG(1), 20, [100, 10, 50, 40, 20, 0], 0, 0.4, 0.5, 0.5)
S2 = runsimulation(StableRNG(1), 20, [100, 10, 50, 40, 20, 0], 0.5, 0.4, 0.5, 0.5)
SC1 = simulationcases(StableRNG(1), 20, [100, 10, 50, 40, 20, 0], 0, 0.4, 0.5, 0.5)
SC2 = simulationcases(StableRNG(1), 20, [100, 10, 50, 40, 20, 0], 0.5, 0.4, 0.5, 0.5)
sim3beta(x) = max(0.1, 2 - 0.1 * x)
rng1 = StableRNG(1)
rng2 = StableRNG(1)
rng3 = StableRNG(1)
sim1 = simulationcases(rng1, 20, [100, 10, 80, 40, 20, 0], 0, 0.5, 0.4, 0.5)
sim2 = simulationcases(rng1, 20, [100, 10, 80, 40, 20, 0], 0.5, 0.5, 0.4, 0.5)
sim3 = simulationcases(rng1, 20, [100, 10, 51, 40, 20, 0], sim3beta, 0.5, 0.4, 0.5)
simtuple1 = ([100, 10, 80, 40, 20, 0], 0, 0.5, 0.4, 0.5, 10)
simtuple2 = ([100, 10, 80, 40, 20, 0], 0.5, 0.5, 0.4, 0.5, 12)
simtuple3 = ([100, 10, 51, 40, 20, 0], sim3beta, 0.5, 0.4, 0.5, nothing)
dict1 = packsimulations(rng2, 20, simtuple1, simtuple2, simtuple3)
gensimtuple1 = packsimulationtuple( ; 
    u0=[100, 10, 80, 40, 20, 0],
    beta=0, 
    eta=0.4, 
    phi=0.5, 
    sigma=0.5, 
    intervention=10,
)
gensimtuple2 = packsimulationtuple( ; 
    u0=[100, 10, 80, 40, 20, 0],
    beta=0.5, 
    eta=0.4, 
    phi=0.5, 
    sigma=0.5, 
    intervention=12,
)
gensimtuple3 = packsimulationtuple( ; 
    u0=[100, 10, 51, 40, 20, 0],
    beta=sim3beta, 
    eta=0.4, 
    phi=0.5, 
    sigma=0.5, 
    intervention=nothing,
)
dict2 = packsimulations(rng3, 20, gensimtuple1, gensimtuple2, gensimtuple3)
expectedS1 = [
    100  10  50  40   20   0
    100   7  21  32   60  11
    100   5   9  26   80  18
    100   3   5  18   94  20
    100   1   3   9  107  22
    100   1   0   7  112  24
    100   0   0   7  113  25
    100   0   0   3  117  25
    100   0   0   3  117  25
    100   0   0   1  119  25
    100   0   0   1  119  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
    100   0   0   0  120  25
]
expectedS2 = [
    100  10  50  40   20   0
     87  20  15  41   57  20
     78  20   8  29   85  27
     72  20   6  21  101  31
     69  14   9  16  112  33
     65  15   5  11  124  35
     62  15   4   7  132  36
     59  11   7   5  138  40
     56   9   7   3  145  41
     55   8   5   3  149  43
     54   8   2   4  152  46
     54   5   1   4  156  47
     53   4   3   3  157  47
     53   2   2   1  162  48
     52   2   1   1  164  48
     52   2   0   2  164  49
     51   2   1   2  164  49
     51   1   2   1  165  49
     51   0   1   2  166  50
     51   0   0   3  166  51
     51   0   0   2  167  51
]
expectedSC1 = [0, 11, 7, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

@testset "simulations u0" begin
    @test simulationu0() == zeros(6)
    @test simulationu0(; S=1) == [1; zeros(Int, 5)]
    @test simulationu0(; N=1) == [zeros(Int, 4); 1; 0]
    @test simulationu0(; S=1, N=1) == [1; zeros(Int, 5)]
    @test_throws ArgumentError simulationu0(; S=2, N=1)
    @test_throws ArgumentError simulationu0(; S=-1)
    @test_throws TypeError simulationu0(; S=1.5)
    @test_throws TypeError simulationu0(; S=1.5, N=2)
    @test simulationu0(; E=1) == [0; 1; zeros(Int, 4)]
    @test simulationu0(; S=1, E=1, N=6) == [1, 1, 0, 0, 4, 0]
    @test_throws ArgumentError simulationu0(; S=1, E=2, N=2)
    @test_throws ArgumentError simulationu0(; E=-1)
    @test_throws TypeError simulationu0(; E=1.5)
    @test simulationu0(; S=1, E=1, I=1, N=6) == [1, 1, 1, 0, 3, 0]
    @test_throws ArgumentError simulationu0(; S=1, E=1, I=2, N=3)
    @test_throws ArgumentError simulationu0(; I=-1)
    @test_throws TypeError simulationu0(; I=1.5)
    @test simulationu0(; S=1, E=1, I=1, Iprime=1) == [1, 1, 1, 1, 0, 0]
    @test simulationu0(; S=1, E=1, I=1, Iprime=1, N=6) == [1, 1, 1, 1, 2, 0]
    @test_throws ArgumentError simulationu0(; S=1, E=1, I=1, Iprime=2, N=4)
    @test_throws ArgumentError simulationu0(; Iprime=-1)
    @test_throws TypeError simulationu0(; Iprime=1.5)
    @test simulationu0(; S=1, E=1, I=1, Iprime=1) == [1, 1, 1, 1, 0, 0]
    @test simulationu0(; S=1, E=1, I=1, Iprime=1, N=5) == [1, 1, 1, 1, 1, 0]
    @test_throws ArgumentError simulationu0(; S=1, E=1, I=1, Iprime=2, N=4)
    @test simulationu0(; R=1) == [zeros(Int, 4); 1; 0]
    @test simulationu0(; R=1, N=1) == [zeros(Int, 4); 1; 0]
    @test_throws ArgumentError simulationu0(; R=1, N=2)
    @test simulationu0(; S=1, E=1, I=1, Iprime=1, R=1) == [1, 1, 1, 1, 1, 0]
    @test simulationu0(; S=1, E=1, I=1, Iprime=1, R=1, N=5) == [1, 1, 1, 1, 1, 0]
    @test_throws ArgumentError simulationu0(; R=-1)
    @test_throws TypeError simulationu0(; R=1.5)
    @test_throws TypeError simulationu0(; N=1.5)
end

@testset "parameters" begin
    @testset "return values with number inputs" begin
        @test RenewalDiD._parameter(0, 1) == 0
        @test RenewalDiD._parameter(0.5, 1) == 0.5
        @test RenewalDiD._parameter(0.5, 2) == 0.5
    end
    @testset "return values with function inputs" begin
        @test RenewalDiD._parameter(x -> x == 1 ? 1 : 0, 1) == 1
        @test RenewalDiD._parameter(x -> x == 1 ? 1 : 0, 2) == 0
    end
    @testset "ensure that `_parameter` always returns a number" begin
        n = RenewalDiD._parameter(x -> sin, 1)
        @test n isa Number
    end
    @testset "throw error for negative parameters" begin
        @test_throws ErrorException RenewalDiD._parameter(-1, 1)
        @test_throws ErrorException RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 1)
        @test_throws ErrorException RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 1, :β)
    end
    @testset "allow possible negative parameters so long as they are never called" begin
        @test RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 2) == 0
        @test RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 2, :β) == 0
    end
    @testset "enforce upper limit on parameter values" begin
        @test RenewalDiD._parameter(x -> 2 * x, 0; upper=1.5) == 0 
        @test RenewalDiD._parameter(x -> 2 * x, 0.5; upper=1.5) == 1 
        @test_throws ErrorException RenewalDiD._parameter(x -> 2 * x, 1; upper=1.5)
    end    
end

@testset "force of infection" begin
    @test RenewalDiD._foi(0, 1, 2, 10, 1) == 0
    @test RenewalDiD._foi(1, 1, 0, 10, 1) == 0.1
    @test RenewalDiD._foi(1, 1, 1, 10, 1) == 0.2
    @test RenewalDiD._foi(1, 1, 2, 10, 1) == 0.3
    @test RenewalDiD._foi(1, 1, 1, 20, 1) == 0.1
    @test RenewalDiD._foi(1, 1, 2, 10, 2) == 0.3
    @test RenewalDiD._foi(x -> x == 2 ? 1 : 0, 1, 2, 10, 1) == 0
    @test RenewalDiD._foi(x -> x == 2 ? 1 : 0, 1, 2, 10, 2) == 0.3
end

@testset "event rates" begin
    @testset "infections" begin
        @test RenewalDiD._simulatedinfections(0, 50, 1, 2, 100, 1) == 0
        @test RenewalDiD._simulatedinfections(1, 50, 1, 2, 100, 1) == 1.5 
    end
    @testset "disease progression" begin
        @test RenewalDiD._diseaseprogression(0, 0, 1) == 0
        @test RenewalDiD._diseaseprogression(10_000, 0.1, 1) == 1000
        @test RenewalDiD._diseaseprogression(10_000, 0, 1) == 0
        @test RenewalDiD._diseaseprogression(10_000, x -> x == 1 ? 1 : 0, 1) == 10_000
        @test RenewalDiD._diseaseprogression(10_000, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._diseaseprogression(10_000, x -> x == 2 ? 1.5 : 0, 1) == 0
    end
    @testset "diagnosis" begin
        @test RenewalDiD._diagnosis(0, 0, 0, 1) == 0
        @test RenewalDiD._diagnosis(500, 0, 0.1, 1) == 0
        @test RenewalDiD._diagnosis(500, 0.9, 0.1, 1) == 50
        @test RenewalDiD._diagnosis(500, 1, x -> x == 1 ? 0.5 : 0, 1) == 500
        @test RenewalDiD._diagnosis(500, 1, x -> x == 1 ? 0.5 : 0, 2) == 0
        @test RenewalDiD._diagnosis(100, 0.5, 0.5, 1) == 50
    end
    @testset "recovery" begin
        @test RenewalDiD._recovery(0, 0, 1) == 0
        @test RenewalDiD._recovery(500, 0.1, 1) == 50
        @test RenewalDiD._recovery(500, x -> x == 1 ? 1 : 0, 1) == 500
        @test RenewalDiD._recovery(500, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._recovery(100, 0.5, 1) == 50
        @test_throws MethodError RenewalDiD._recovery(1, 1, 0, 1)
    end
    @testset "calculate population size" begin
        @test RenewalDiD._n_seir(zeros(Int, 6)) == 0
        @test RenewalDiD._n_seir(ones(Int, 6)) == 5
        @test RenewalDiD._n_seir([3, 2, 1, 0, 0, 0]) == 6
        @test RenewalDiD._n_seir([3, 2, 1, 0, 0, 4]) == 6
        @test_throws ArgumentError RenewalDiD._n_seir(ones(Int, 7))
        @test_throws MethodError RenewalDiD._n_seir(ones(6))
        @test_throws ArgumentError RenewalDiD._n_seir([1, 1, 0, 0, -1, 0])
    end
    @testset "vector of event rates" begin       
        u0_1 = simulationu0(; N=1)
        @test testseirrates(u0_1, 1) == zeros(5)
        # vector of length 5 for u0
        @test_throws ArgumentError testseirrates([0, 0, 0, 1, 0], 1) 
        u0_2 = simulationu0(; S=10_000, I=500, N=20_000)
        @test testseirrates(u0_2, 1; beta=0.5) == [125, 0, 0, 0, 0]
        u0_3 = simulationu0(; S=10_000, I=500, Iprime=500, N=20_000)
        @test testseirrates(u0_3, 1; beta=0.5) == [250, 0, 0, 0, 0]
        u0_4 = simulationu0(; S=10_000, I=1000, Iprime=500, N=20_000)
        @test testseirrates(u0_4, 1; beta=0.5) == [375, 0, 0, 0, 0]
        u0_5 = simulationu0(; S=10_000, E=1000, I=500, Iprime=500, N=20_000)
        @test testseirrates(u0_5, 1; beta=0.5) == [250, 0, 0, 0, 0]
        @test testseirrates(u0_5, 1; beta=0.5, eta=0, phi=0.5) == [250, 0, 0, 0, 0]
        @test testseirrates(u0_5, 1; beta=0.5, eta=1, phi=0.5) == [250, 0, 500, 500, 500]
        u0_8 = simulationu0(; S=10_000, E=1000, I=1000, N=20_000)
        r8 = testseirrates(u0_8, 1; beta=0.5, sigma=0.5, eta=0.6, phi=0.5)
        @test r8 == [250, 500, 600, 600, 0]
        u0_9 = simulationu0(; S=10_000, E=1000, I=600, Iprime=400, N=20_000)
        r9 = testseirrates(u0_9, 1; beta=0.5, sigma=0.5, eta=0.6, phi=0.5)
        @test r9 == [250, 500, 360, 360, 240]
        u0_11 = simulationu0(; S=10_000, E=1000, I=600, Iprime=1000, N=20_000)
        r11 = testseirrates(u0_11, 1; beta=0.5, sigma=0.5, eta=0.6, phi=0.5)
        @test r11 == [400, 500, 360, 360, 600]
    end
end

@testset "next event" begin
    @testset "time step" begin
        @test RenewalDiD._tstep(StableRNG(1), [80, 10, 50, 30, 20]) == expectedv1
    end
    @testset "identifying next event" begin
        @test RenewalDiD._nextevent(Random.default_rng(), [100, 0, 0, 0, 0]) == 1
        @test RenewalDiD._nextevent(StableRNG(1), [0, 0, 100, 0, 0]) == 3
        @test RenewalDiD._nextevent(StableRNG(1), [100, 10, 50, 30, 60]) == 3
        @test RenewalDiD._nextevent(StableRNG(2), [100, 10, 50, 30, 60]) == 5
    end
    @testset "update event" begin
        # this test set mutates the same vector multiple times so if multiple tests fail it 
        # may only be a problem with the first of the failing functions
        u1 = [100, 100, 100, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 1)  # infection
        @test u1 == [99, 101, 100, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 2)  # disease progression,
        @test u1 == [99, 100, 101, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 3)  # diagnosis
        @test u1 == [99, 100, 100, 101, 100, 1]
        RenewalDiD._updateevent!(u1, 4)  # recovery from I
        @test u1 == [99, 100, 99, 101, 101, 1]
        RenewalDiD._updateevent!(u1, 5)  # recovery from Iprime
        @test u1 == [99, 100, 99, 100, 102, 1]
    end
end

@testset "simulate a day" begin
    # some functions in this test set mutate the same vector multiple times so if multiple 
    # tests fail it may only be a problem with the first of the failing functions
    u2 = simulationu0(; S=1, Iprime=1)
    testsimulateday!(StableRNG(1), u2, 1)
    @test u2 == [1, 0, 0, 1, 0, 0]
    testsimulateday!(StableRNG(1), u2, 1; beta=100)
    @test u2 == [0, 1, 0, 1, 0, 0]
    u3 = simulationu0(; E=1)
    testsimulateday!(StableRNG(1), u3, 1; sigma=100)
    @test u3 == [0, 0, 1, 0, 0, 0]
    testsimulateday!(StableRNG(1), u3, 1; eta=1000)
    @test u3 == [0, 0, 0, 0, 1, 0]
    u4 = simulationu0(; E=1)
    testsimulateday!(StableRNG(1), u4, 1; sigma=100)
    @test u4 == [0, 0, 1, 0, 0, 0]
    testsimulateday!(StableRNG(1), u4, 2; eta=1000, phi=0.99999)
    @test u4 == [0, 0, 0, 0, 1, 1]  # diagnosed then recovered
    testsimulateday!(StableRNG(1), u4, 3; eta=1000)
    @test u4 == [0, 0, 0, 0, 1, 1]
    @test_throws ErrorException testsimulateday!(StableRNG(1), u4, 2; eta=1000, phi=10)
end

@testset "run whole simulation" begin
    @test_throws DimensionMismatch runsimulation(
        20, [100, 10, 50, 30, 40], 0.5, 0.4, 0.5, 0.4
    )
    u0 = simulationu0(; S=100, E=10, I=50, Iprime=40, R=20)
    @test_throws MethodError runsimulation(20.2, u0, 0.5, 0.4, 0.5, 0.5)
    @test_throws ErrorException runsimulation(20, u0, x -> 0.5 - 0.1 * x, 0.4, 0.5, 0.5)
    @test_throws ArgumentError runsimulation(-20, u0, 0, 0.4, 0.5, 0.5)
    @test S1 == expectedS1
    @test S2 == expectedS2
    @test simulationcases(zeros(10, 6)) == zeros(10)
    @test simulationcases(zeros(12, 6)) == zeros(12)
    @test_throws ArgumentError simulationcases(zeros(12, 5))
    @test simulationcases(S1) == expectedSC1
    @test SC1 == expectedSC1
    @test SC2 isa Vector{Int}
    @test length(SC2) == 21
end 

@testset "group simulations" begin
    @test dict1.interventions == InterventionMatrix(20, 10, 12, nothing)
    @test dict1.Ns == [250, 250, 221]
    @testset for i in 1:3 
        @test dict1.observedcases[:, i] == [sim1, sim2, sim3][i]
    end
    @test gensimtuple1 == simtuple1
    @test gensimtuple2 == simtuple2
    @test gensimtuple3 == simtuple3
    @test dict1 == dict2
end
