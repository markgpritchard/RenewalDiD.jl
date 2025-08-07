# test functions that generate simulations 

using RenewalDiD
using StableRNGs
using Test

expectedv1 = 0.002143243061437003
S1 = runsimulation(
    StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5
)
S2 = runsimulation(
    StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5
)
SC1 = simulationcases(
    StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5
)
SC2 = simulationcases(
    StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5
)
sim3beta(x) = max(0.1, 2 - 0.1 * x)
rng1 = StableRNG(1)
rng2 = StableRNG(1)
rng3 = StableRNG(1)
sim1 = simulationcases(
    rng1, 20, [100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5
)
sim2 = simulationcases(
    rng1, 20, [100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5
)
sim3 = simulationcases(
    rng1, 20, [100, 10, 50, 1, 40, 20, 0], sim3beta, 0.4, 0.5, 0.5, 0.5
)
simtuple1 = ([100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5, 10)
simtuple2 = ([100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5, 12)
simtuple3 = ([100, 10, 50, 1, 40, 20, 0], sim3beta, 0.4, 0.5, 0.5, 0.5, nothing)
dict1 = packsimulations(rng2, 20, simtuple1, simtuple2, simtuple3)
gensimtuple1 = packsimulationtuple( ; 
    u0=[100, 10, 50, 30, 40, 20, 0],
    beta=0, 
    gamma=0.4, 
    delta=0.5, 
    theta=0.5, 
    sigma=0.5, 
    intervention=10,
)
gensimtuple2 = packsimulationtuple( ; 
    u0=[100, 10, 50, 30, 40, 20, 0],
    beta=0.5, 
    gamma=0.4, 
    delta=0.5, 
    theta=0.5, 
    sigma=0.5, 
    intervention=12,
)
gensimtuple3 = packsimulationtuple( ; 
    u0=[100, 10, 50, 1, 40, 20, 0],
    beta=sim3beta, 
    gamma=0.4, 
    delta=0.5, 
    theta=0.5, 
    sigma=0.5, 
    intervention=nothing,
)
dict2 = packsimulations(rng3, 20, gensimtuple1, gensimtuple2, gensimtuple3)
expectedS1 = [
    100  10  50  30  40   20   0
    100   6  37  20   9   78  11
    100   4  25  13   3  105  18
    100   2  15   9   1  123  24
    100   0   8   9   0  133  26
    100   0   7   7   2  134  28
    100   0   6   7   0  137  28
    100   0   3   5   1  141  30
    100   0   2   4   2  142  31
    100   0   0   4   0  146  31
    100   0   0   3   0  147  32
    100   0   0   3   0  147  32
    100   0   0   1   0  149  34
    100   0   0   1   0  149  34
    100   0   0   1   0  149  34
    100   0   0   1   0  149  34
    100   0   0   1   0  149  34
    100   0   0   0   0  150  35
    100   0   0   0   0  150  35
    100   0   0   0   0  150  35
    100   0   0   0   0  150  35
]
expectedS2 = [
    100  10  50  30  40   20   0
     87  16  35  20  14   78  12
     76  20  22  11   8  113  24
     73  16  18  10   4  129  29
     69  18  10   8   1  144  33
     64  14  12   6   1  153  39
     61  11  10   8   1  159  41
     60   8  11   6   0  165  44
     58   3   8  10   2  169  46
     57   3   6   6   1  177  50
     56   3   5   5   0  181  52
     52   6   3   5   1  183  53
     52   4   4   5   1  184  54
     52   3   3   3   2  187  57
     50   3   2   3   1  191  59
     49   4   2   3   0  192  59
     49   3   3   1   0  194  61
     49   2   2   0   0  197  62
     49   1   2   1   0  197  62
     49   0   3   1   0  197  62
     49   0   3   0   0  198  63
]
expectedSC1 = [0, 11, 7, 6, 2, 2, 0, 2, 1, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0]

@testset "simulations u0" begin
    @test simulationu0() == zeros(7)
    @test simulationu0(; s=1) == [1; zeros(Int, 6)]
    @test simulationu0(; n=1) == [zeros(Int, 5); 1; 0]
    @test simulationu0(; s=1, n=1) == [1; zeros(Int, 6)]
    @test_throws ArgumentError simulationu0(; s=2, n=1)
    @test_throws ArgumentError simulationu0(; s=-1)
    @test_throws TypeError simulationu0(; s=1.5)
    @test_throws TypeError simulationu0(; s=1.5, n=2)
    @test simulationu0(; e=1) == [0; 1; zeros(Int, 5)]
    @test simulationu0(; s=1, e=1, n=6) == [1, 1, 0, 0, 0, 4, 0]
    @test_throws ArgumentError simulationu0(; s=1, e=2, n=2)
    @test_throws ArgumentError simulationu0(; e=-1)
    @test_throws TypeError simulationu0(; e=1.5)
    @test simulationu0(; s=1, e=1, i_n=1, n=6) == [1, 1, 1, 0, 0, 3, 0]
    @test_throws ArgumentError simulationu0(; s=1, e=1, i_n=2, n=3)
    @test_throws ArgumentError simulationu0(; i_n=-1)
    @test_throws TypeError simulationu0(; i_n=1.5)
    @test simulationu0(; s=1, e=1, i_n=1, i_f=1) == [1, 1, 1, 1, 0, 0, 0]
    @test simulationu0(; s=1, e=1, i_n=1, i_f=1, n=6) == [1, 1, 1, 1, 0, 2, 0]
    @test_throws ArgumentError simulationu0(; s=1, e=1, i_n=1, i_f=2, n=4)
    @test_throws ArgumentError simulationu0(; i_f=-1)
    @test_throws TypeError simulationu0(; i_f=1.5)
    @test simulationu0(; s=1, e=1, i_n=1, i_f=1, i_d=1) == [1, 1, 1, 1, 1, 0, 0]
    @test simulationu0(; s=1, e=1, i_n=1, i_f=1, i_d=1, n=6) == [1, 1, 1, 1, 1, 1, 0]
    @test_throws ArgumentError simulationu0(; s=1, e=1, i_n=1, i_f=1, i_d=2, n=5)
    @test_throws ArgumentError simulationu0(; i_d=-1)
    @test_throws TypeError simulationu0(; i_d=1.5)
    @test simulationu0(; r=1) == [zeros(Int, 5); 1; 0]
    @test simulationu0(; r=1, n=1) == [zeros(Int, 5); 1; 0]
    @test_throws ArgumentError simulationu0(; r=1, n=2)
    @test simulationu0(; s=1, e=1, i_n=1, i_f=1, i_d=1, r=1) == [1, 1, 1, 1, 1, 1, 0]
    @test simulationu0(; s=1, e=1, i_n=1, i_f=1, i_d=1, r=1, n=6) == [1, 1, 1, 1, 1, 1, 0]
    @test_throws ArgumentError simulationu0(; r=-1)
    @test_throws TypeError simulationu0(; r=1.5)
    @test_throws TypeError simulationu0(; n=1.5)
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
        @test RenewalDiD._parameter(x -> 2 * x, 0; upper=1) == 0 
        @test RenewalDiD._parameter(x -> 2 * x, 0.5; upper=1) == 1 
        @test_throws ErrorException RenewalDiD._parameter(x -> 2 * x, 1; upper=1)
    end    
end

@testset "force of infection" begin
    @test RenewalDiD._foi(0, 1, 1, 1, 10, 1) == 0
    @test RenewalDiD._foi(1, 1, 0, 0, 10, 1) == 0.1
    @test RenewalDiD._foi(1, 1, 1, 0, 10, 1) == 0.2
    @test RenewalDiD._foi(1, 1, 1, 1, 10, 1) == 0.3
    @test RenewalDiD._foi(1, 1, 1, 0, 20, 1) == 0.1
    @test RenewalDiD._foi(1, 1, 1, 1, 10, 2) == 0.3
    @test RenewalDiD._foi(x -> x == 2 ? 1 : 0, 1, 1, 1, 10, 1) == 0
    @test RenewalDiD._foi(x -> x == 2 ? 1 : 0, 1, 1, 1, 10, 2) == 0.3
end

@testset "event rates" begin
    @testset "infections" begin
        @test RenewalDiD._simulatedinfections(0, 50, 1, 1, 1, 100, 1) == 0
        @test RenewalDiD._simulatedinfections(1, 50, 1, 1, 1, 100, 1) == 1.5 
    end
    @testset "disease progression (never diagnosed)" begin
        @test RenewalDiD._diseaseprogression_to_i_n(0, 0, 0, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0, 0.1, 1) == 1000
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 1, 0.1, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0.25, 0.1, 1) == 750
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_n(1000, 1.5, 0.1, 1)
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, x -> x == 1 ? 1 : 0, 1, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, x -> x == 1 ? 1 : 0, 1, 2) == 10_000
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, x -> x == 2 ? 1.5 : 0, 1, 1) == 10_000
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_n(
            10_000, x -> x == 2 ? 1.5 : 0, 1, 2
        )
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0, x -> x == 1 ? 1 : 0, 1) == 10_000
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0, x -> x == 1 ? 1 : 0, 2) == 0
    end
    @testset "disease progression (will be diagnosed)" begin
        @test RenewalDiD._diseaseprogression_to_i_f(0, 0, 0, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 0, 0.1, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 1, 0.1, 1) == 1000
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 0.25, 0.1, 1) == 250
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_f(1000, 1.5, 0.1, 1)
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, x -> x == 1 ? 1 : 0, 1, 1) == 10_000
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, x -> x == 1 ? 1 : 0, 1, 2) == 0
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, x -> x == 2 ? 1.5 : 0, 1, 1) == 0
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_f(
            10_000, x -> x == 2 ? 1.5 : 0, 1, 2
        )
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 1, x -> x == 1 ? 1 : 0, 1) == 10_000
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 1, x -> x == 1 ? 1 : 0, 2) == 0
    end
    @testset "diagnosis" begin
        @test RenewalDiD._diagnosis(0, 0, 1) == 0
        @test RenewalDiD._diagnosis(500, 0.1, 1) == 50
        @test RenewalDiD._diagnosis(500, x -> x == 1 ? 1 : 0, 1) == 500
        @test RenewalDiD._diagnosis(500, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._diagnosis(100, 0.5, 1) == 50
    end
    @testset "recovery" begin
        @test RenewalDiD._recovery(0, 0, 1) == 0
        @test RenewalDiD._recovery(500, 0.1, 1) == 50
        @test RenewalDiD._recovery(500, x -> x == 1 ? 1 : 0, 1) == 500
        @test RenewalDiD._recovery(500, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._recovery(100, 0.5, 1) == 50
        @test_throws ErrorException RenewalDiD._recovery(1, 1, 0, 1)
        @test RenewalDiD._recovery(0, 0.1, 0.5, 1) == 0
        @test RenewalDiD._recovery(500, 0.1, Inf, 1) == 50
        @test RenewalDiD._recovery(500, 0.1, 0.2, 1) == 100
        @test RenewalDiD._recovery(500, x -> 0.1 * x, 0.2, 0) == 0
        @test RenewalDiD._recovery(500, x -> 0.1 * x, 0.2, 1) == 100
        @test_throws ErrorException RenewalDiD._recovery(500, x -> 0.1 * x, 0.2, 2)
        @test RenewalDiD._recovery(500, x -> 0.1 * x, x -> x == 3 ? 0.2 : Inf, 0) == 0
        @test RenewalDiD._recovery(500, x -> 0.1 * x, x -> x == 3 ? 0.2 : Inf, 1) == 50
        @test RenewalDiD._recovery(500, x -> 0.1 * x, x -> x == 3 ? 0.2 : Inf, 2) == 100
        @test_throws ErrorException RenewalDiD._recovery(
            500, x -> 0.1 * x, x -> x == 3 ? 0.2 : Inf, 3
        )
    end
    @testset "calculate population size" begin
        @test RenewalDiD._n_seir(zeros(Int, 7)) == 0
        @test RenewalDiD._n_seir(ones(Int, 7)) == 6
        @test_throws ArgumentError RenewalDiD._n_seir(ones(Int, 8))
        @test_throws MethodError RenewalDiD._n_seir(ones(7))
        @test_throws ArgumentError RenewalDiD._n_seir([1, 1, 0, 0, 0, -1, 0])
    end
    @testset "vector of event rates" begin
        u0_1 = simulationu0(; n=1)
        r1 = RenewalDiD._seirrates(u0_1, 1, 0, 0, 0, 0, 0)
        @test r1 == zeros(6)
        @test_throws ArgumentError RenewalDiD._seirrates(  # vector of length 6 for u0
            [0, 0, 0, 0, 1, 0], 1, 0, 0, 0, 0, 0
        ) 
        u0_2 = simulationu0(; s=10_000, i_n=500, n=20_000)
        r2 = RenewalDiD._seirrates(u0_2, 1, 0.5, 0, 0, 0, 0)
        @test r2 == [125, 0, 0, 0, 0, 0]
        u0_3 = simulationu0(; s=10_000, i_n=500, i_f=500, n=20_000)
        r3 = RenewalDiD._seirrates(u0_3, 1, 0.5, 0, 0, 0, 0)
        @test r3 == [250, 0, 0, 0, 0, 0]
        u0_4 = simulationu0(; s=10_000, i_n=500, i_f=500, i_d=500, n=20_000)
        r4 = RenewalDiD._seirrates(u0_4, 1, 0.5, 0, 0, 0, 0)
        @test r4 == [375, 0, 0, 0, 0, 0]
        u0_5 = simulationu0(; s=10_000, e=1000, i_n=500, i_f=500, n=20_000)
        r5 = RenewalDiD._seirrates(u0_5, 1, 0.5, 0, 0, 0, 0)
        @test r5 == [250, 0, 0, 0, 0, 0]
        r6 = RenewalDiD._seirrates(u0_5, 1, 0.5, 0, 0, 0, 0.5)
        @test r6 == [250, 500, 0, 0, 0, 0]
        r7 = RenewalDiD._seirrates(u0_5, 1, 0.5, 0, 0, 1, 0.5)
        @test r7 == [250, 0, 500, 0, 0, 0]
        u0_8 = simulationu0(; s=10_000, e=1000, i_n=1000, n=20_000)
        r8 = RenewalDiD._seirrates(u0_8, 1, 0.5, 0, 0.5, 0.6, 0.5)
        @test r8 == [250, 200, 300, 0, 0, 0]
        u0_9 = simulationu0(; s=10_000, e=1000, i_n=600, i_f=400, n=20_000)
        r9 = RenewalDiD._seirrates(u0_9, 1, 0.5, 0, 0.5, 0.6, 0.5)
        @test r9 == [250, 200, 300, 200, 0, 0]
        r10 = RenewalDiD._seirrates(u0_9, 1, 0.5, 1/3, 0.5, 0.6, 0.5)
        @test r10 == [250, 200, 300, 200, 200, 0]
        u0_11 = simulationu0(; s=10_000, e=1000, i_n=600, i_f=400, i_d=1000, n=20_000)
        r11 = RenewalDiD._seirrates(u0_11, 1, 0.5, 1/3, 0.5, 0.6, 0.5)
        @test r11 == [500, 200, 300, 200, 200, 1000]
    end
end

@testset "next event" begin
    @testset "time step" begin
        @test RenewalDiD._tstep(StableRNG(1), [80, 10, 50, 30, 60, 20]) == expectedv1
    end
    @testset "next event time" begin  
        # function not used in the package so no longer needs testing (and is deleted)
    end
    @testset "identifying next event" begin
        @test RenewalDiD._nextevent([100, 0, 0, 0, 0, 0]) == 1
        @test RenewalDiD._nextevent(StableRNG(1), [0, 0, 100, 0, 0, 0]) == 3
        @test RenewalDiD._nextevent(StableRNG(1), [100, 10, 50, 30, 60, 0]) == 3
        @test RenewalDiD._nextevent(StableRNG(2), [100, 10, 50, 30, 60, 0]) == 5
    end
    @testset "update event" begin
        # this test set mutates the same vector multiple times so if multiple tests fail it 
        # may only be a problem with the first of the failing functions
        u1 = [100, 100, 100, 100, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 1)  # infection
        @test u1 == [99, 101, 100, 100, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 2)  # disease progression, subset who won't be diagnosed
        @test u1 == [99, 100, 101, 100, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 3)  # disease progression, subset who will be diagnosed
        @test u1 == [99, 99, 101, 101, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 4)  # diagnosis
        @test u1 == [99, 99, 101, 100, 101, 100, 1]
        RenewalDiD._updateevent!(u1, 5)  # recovery from i_n
        @test u1 == [99, 99, 100, 100, 101, 101, 1]
        RenewalDiD._updateevent!(u1, 6)  # recovery from i_d
        @test u1 == [99, 99, 100, 100, 100, 102, 1]
    end
end

@testset "simulate a day" begin
    # some functions in this test set mutate the same vector multiple times so if multiple 
    # tests fail it may only be a problem with the first of the failing functions
    u2 = simulationu0(; s=1, i_f=1)
    RenewalDiD._simulateday!(StableRNG(1), u2, 1, 0, 0, 0, 0, 0)
    @test u2 == [1, 0, 0, 1, 0, 0, 0]
    RenewalDiD._simulateday!(StableRNG(1), u2, 1, 100, 0, 0, 0, 0)
    @test u2 == [0, 1, 0, 1, 0, 0, 0]
    u3 = simulationu0(; e=1)
    [0, 1, 0, 0, 0, 0, 0]
    RenewalDiD._simulateday!(StableRNG(1), u3, 1, 0, 0, 0, 0, 100)
    @test u3 == [0, 0, 1, 0, 0, 0, 0]
    RenewalDiD._simulateday!(StableRNG(1), u3, 1, 0, 100, 1000, 0, 0)
    @test u3 == [0, 0, 0, 0, 0, 1, 0]
    u4 = simulationu0(; e=1)
    RenewalDiD._simulateday!(StableRNG(1), u4, 1, 0, 0, 0, 1, 100)
    @test u4 == [0, 0, 0, 1, 0, 0, 0]
    RenewalDiD._simulateday!(StableRNG(1), u4, 1, 0, 0, 100, 1, 0)
    @test u4 == [0, 0, 0, 0, 1, 0, 1]
    RenewalDiD._simulateday!(StableRNG(1), u4, 1, 0, 100, 1000, 0, 0)
    @test u4 == [0, 0, 0, 0, 0, 1, 1]
end

@testset "run whole simulation" begin
    @test_throws DimensionMismatch runsimulation(
        20, [100, 10, 50, 30, 40, 20], 0.5, 0.4, 0.5, 0.4, 0.5
    )
    u0 = simulationu0(; s=100, e=10, i_n=50, i_d=30, i_f=40, r=20)
    @test_throws MethodError runsimulation(20.2, u0, 0.5, 0.4, 0.5, 0.4, 0.5)
    @test_throws ErrorException runsimulation(20, u0, x -> 0.5 - 0.1 * x, 0.4, 0.5, 0.4, 0.5)
    @test_throws ErrorException runsimulation(20, u0, 1.5, 0.4, 0.5, x -> 0.5 + 0.1 * x, 0.5)
    @test_throws ArgumentError runsimulation(-20, u0, 0, 0.4, 0.5, 0.5, 0.5)
    @test S1 == expectedS1
    @test S2 == expectedS2
    @test simulationcases(zeros(10, 7)) == zeros(10)
    @test simulationcases(zeros(12, 7)) == zeros(12)
    @test_throws ArgumentError simulationcases(zeros(12, 6))
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
