using RenewalDiD
using Test
using StableRNGs

@testset "RenewalDiD.jl" begin

@testset "intervention matrix" begin
    @test InterventionMatrix{Int}(2, 2; mutewarnings=true) isa AbstractMatrix
    # behaviour of `; mutewarnings` is tested below
    @test size(InterventionMatrix{Int}(2, 2; mutewarnings=true)) == (2, 1)
    @test size(InterventionMatrix{Int}(3, 2; mutewarnings=true)) == (3, 1)
    @test size(InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true)) == (3, 2)
    M1 = InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true) 
    @test M1[1] == 0
    @test M1[1, 1] == 0
    @test M1[3, 1] == 1
    @test M1[3] == 1
    @test M1[4] == 0
    @test M1[6] == 1
    @test_throws BoundsError M1[7] 
    M2 = InterventionMatrix{Int}(3, [2, nothing])
    @test M2[1, 2] == 0
    @test M2[3, 2] == 0
    @test InterventionMatrix(2, 2; mutewarnings=true) == InterventionMatrix{Int}(
        2, 2; mutewarnings=true
    )  
    # note that this only tests whether the contructor without a type works, not that it 
    # generates a type `Int`; tests of `Base.show` below demonstrate that the default type
    # is `Int`
    @test InterventionMatrix(3, 2; mutewarnings=true) == InterventionMatrix{Int}(
        3, 2; mutewarnings=true
    )
    @test InterventionMatrix(3, [2, 3]; mutewarnings=true) == M1
    @test RenewalDiD._duration(M1) == 3
    @test RenewalDiD._duration(InterventionMatrix(4, [2, 3]; mutewarnings=true)) == 4
    @test_throws BoundsError M1[4, 1]
    @test_throws BoundsError M1[1, 3]
    @test_throws MethodError InterventionMatrix(3.3, 2)
    @test_throws MethodError InterventionMatrix(3, 2.3)
    @test_throws MethodError InterventionMatrix(3, [2.3, 2])
    @test InterventionMatrix{Bool}(3, [2, 3]; mutewarnings=true)[3, 1]
    @test InterventionMatrix(3, 2, 3; mutewarnings=true) == M1
    @test InterventionMatrix(4, 2, 3; mutewarnings=true) == InterventionMatrix{Int}(
        4, [2, 3]; mutewarnings=true
    ) 
    @test InterventionMatrix(4, 2, 3, 4; mutewarnings=true) == InterventionMatrix{Int}(
        4, [2, 3, 4]; mutewarnings=true
    )
    @test InterventionMatrix(4, 2, 3, nothing) == InterventionMatrix{Int}(4, [2, 3, 5])
    tw1 = "InterventionMatrix with no intervention in any group"
    @test_warn tw1 InterventionMatrix(2, nothing, nothing)
    @test_nowarn InterventionMatrix{Int}(2, [2, 3])
    @test_nowarn InterventionMatrix(4, 2, 3, nothing)
    @test_nowarn InterventionMatrix(2, nothing, nothing; mutewarnings=true)
    @test_warn tw1 InterventionMatrix(2, nothing, nothing; mutewarnings=false)
    tw2 = "All groups in InterventionMatrix have intervention before end of duration"
    @test_warn tw2 InterventionMatrix(3, 2, 3)
    @test_nowarn InterventionMatrix(3, 2, 3; mutewarnings=true)
    @test_warn tw2 InterventionMatrix(3, 2, 3; mutewarnings=false)
    M3 = InterventionMatrix(10, 4, 5, 8, 8, nothing)
    @test RenewalDiD._showtimes(M1) == [1, 2, 3]
    @test RenewalDiD._showtimes(M3) == [1, 4, 5, 8, 10]
    @test RenewalDiD._showcontents(M1, RenewalDiD._showtimes(M1)) == [
        0  0
        1  0 
        1  1
    ]
    @test RenewalDiD._showcontents(M3, RenewalDiD._showtimes(M3)) == [
        0  0  0  0  0
        1  0  0  0  0
        1  1  0  0  0
        1  1  1  1  0
        1  1  1  1  0
    ]
    @testset "show intervention matrix M1" begin
        ts, cs = RenewalDiD._showstrings(M1)
        @test ts == ["1", "2", "3"]
        @test cs == [
            "0"  "0"
            "1"  "0" 
            "1"  "1"
        ]
        @test RenewalDiD._showcombinedstrings(M1) == [
            "1"  "0"  "0"
            "2"  "1"  "0" 
            "3"  "1"  "1"
        ]
        @test RenewalDiD._showheader(M1) == ["time", "1", "2"]
        O1 = "3×2 InterventionMatrix{Int64}\n time │ 1  2\n──────┼──────\n    1 │ 0  0\n    2 │ 1  0\n    3 │ 1  1\n"
        @test repr("text/plain", M1) == O1
    end 
    @testset "show intervention matrix M3" begin
        ts, cs = RenewalDiD._showstrings(M3)
        @test ts == ["1", "⋮", "4", "5", "⋮", "8", "⋮", "10"]
        @test cs == [
            "0"  "0"  "0"  "0"  "0"
            "⋮"  "⋮"  "⋮"  "⋮"  "⋮"
            "1"  "0"  "0"  "0"  "0"
            "1"  "1"  "0"  "0"  "0"
            "⋮"  "⋮"  "⋮"  "⋮"  "⋮"
            "1"  "1"  "1"  "1"  "0"
            "⋮"  "⋮"  "⋮"  "⋮"  "⋮"
            "1"  "1"  "1"  "1"  "0"
        ]
    end
    @testset "show intervention matrix M4" begin
        M4 = InterventionMatrix{Float64}(5, 2, 3, nothing)   
        ts, cs = RenewalDiD._showstrings(M4)
        @test cs == [
            "0.0"  "0.0"  "0.0" 
            "1.0"  "0.0"  "0.0" 
            "1.0"  "1.0"  "0.0"
            "⋮"  "⋮"  "⋮"
            "1.0"  "1.0"  "0.0"
        ]
        @test RenewalDiD._showcombinedstrings(M4) == [
            "1"  "0.0"  "0.0"  "0.0" 
            "2"  "1.0"  "0.0"  "0.0" 
            "3"  "1.0"  "1.0"  "0.0"
            "⋮"  "⋮"  "⋮"  "⋮"
            "5"  "1.0"  "1.0"  "0.0"
        ]
        @test RenewalDiD._showheader(M4) == ["time", "1", "2", "3"]
        O4 = "5×3 InterventionMatrix{Float64}\n time │   1    2    3\n──────┼───────────────\n    1 │ 0.0  0.0  0.0\n    2 │ 1.0  0.0  0.0\n    3 │ 1.0  1.0  0.0\n    ⋮ │   ⋮    ⋮    ⋮\n    5 │ 1.0  1.0  0.0\n"
        @test repr("text/plain", M4) == O4
    end
end  # @testset "intervention matrix"

@testset "generation interval" begin
    @test g_covid(0) == 0
    @test g_covid(1) == 0.0440204506
    @test_throws MethodError g_covid(1.2)
    @test g_covid(30) == 0.0002057625
    @test g_covid(31) == 0
    @test g_covid(-1) == 0
    @test g_seir(0; gamma=0.4, sigma=0.5) == 0
    @test g_seir(1; gamma=0.4, sigma=0.5) == 0.12757877264601186
    @test g_seir(1; gamma=0.5, sigma=0.5) == 0.15163266492815836
    @test g_seir(1; gamma=0.5) == 0.15163266492815836
    @test g_seir(0; gamma=0.5) == 0
    @test_throws UndefKeywordError g_seir(1)
    @test g_seir(-1; gamma=0.5) == 0
    @test g_seir(-1; gamma=0.4, sigma=0.5) == 0
    v1 = [g_covid(x) for x in -1000:1:1000]
    @test sum(v1) <= 1
    @test minimum(v1) >= 0
    v2 = [g_seir(x; gamma=0.5) for x in -1000:1:1000]
    @test sum(v2) <= 1
    @test minimum(v2) >= 0
    v3 = [g_seir(x; gamma=0.4, sigma=0.5) for x in -1000:1:1000]
    @test sum(v3) <= 1
    @test minimum(v3) >= 0
    @test generationtime([0.0, 0.2, 0.4], 0) == 0
    @test generationtime([0.0, 0.2, 0.4], 1) == 0.2
    @test generationtime([0.0, 0.2, 0.4], -1) == 0
    @test generationtime([0.0, 0.2, 0.4], 3) == 0
    @test_throws MethodError generationtime([0.0, 0.2, 0.4], 2.3)
    @test_throws ArgumentError generationtime([0.0, 0.2, 0.4, 0.8], 0)
    @test_throws ArgumentError generationtime([0.0, 0.2, 0.4, -0.1], 0)
    @test generationtime(x -> x < 2 ? 0.0 : 1/x^2, 0) == 0
    @test generationtime(x -> x < 2 ? 0.0 : 1/x^2, 2) == 0.25
    @test generationtime(x -> x == -1 ? 0.25 : 0, -1) == 0
    @test_throws MethodError generationtime(x -> x < 2 ? 0.0 : 1/x^2, 1.5)
    @test_throws ArgumentError generationtime(x -> 1/x^2, 0)
    @test_throws ArgumentError generationtime(x -> x == 12 ? -0.5 : 0.0, 0)
    @test generationtime(x -> x == 12 ? -0.5 : 0.0, 1; t_max=10) == 0
    @test generationtime(x -> x == 12 ? -0.5 : 0.0, 12; t_max=10) == 0
    @testset "vector g_seir" begin
        v1 = vectorg_seir(0.4, 0.5)
        @test generationtime(v1, 1) == g_seir(1; gamma=0.4, sigma=0.5)
        @test generationtime(v1, 12) == g_seir(12; gamma=0.4, sigma=0.5)
        @test generationtime(v1, 30) == 0
        v2 = vectorg_seir(0.4, 0.5; t_max=10)
        @test generationtime(v2, 2) == g_seir(2; gamma=0.4, sigma=0.5)
        @test generationtime(v2, 12) == 0
        v3 = vectorg_seir(0.4)
        @test generationtime(v3, 1) == g_seir(1; gamma=0.4)
        @test generationtime(v3, 12) == g_seir(12; gamma=0.4, sigma=0.4)
        @test generationtime(v3, 30) == 0
        v4 = vectorg_seir(0.4; t_max=10)
        @test generationtime(v4, 2) == g_seir(2; gamma=0.4)
        @test generationtime(v4, 12) == 0
    end
end  # @testset "generation interval"

@testset "simulations" begin
    @testset "parameters" begin
        @test RenewalDiD._parameter(0, 1) == 0
        @test RenewalDiD._parameter(0.5, 1) == 0.5
        @test RenewalDiD._parameter(0.5, 2) == 0.5
        @test RenewalDiD._parameter(x -> x == 1 ? 1 : 0, 1) == 1
        @test RenewalDiD._parameter(x -> x == 1 ? 1 : 0, 2) == 0
        @test_throws ErrorException RenewalDiD._parameter(-1, 1)
        @test_throws ErrorException RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 1)
        @test RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 2) == 0
        @test_throws ErrorException RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 1, :β)
        @test RenewalDiD._parameter(x -> x == 1 ? -1 : 0, 2, :β) == 0
        @test RenewalDiD._parameter(nothing, 2, :β) == 0 
    end
    @testset "force of infection" begin
        @test RenewalDiD._foi(0, 1, 1, 10, 1) == 0
        @test RenewalDiD._foi(1, 1, 0, 10, 1) == 0.1
        @test RenewalDiD._foi(1, 1, 1, 10, 1) == 0.2
        @test RenewalDiD._foi(1, 1, 1, 20, 1) == 0.1
        @test RenewalDiD._foi(1, 1, 1, 10, 2) == 0.2
        @test RenewalDiD._foi(x -> x == 2 ? 1 : 0, 1, 1, 10, 1) == 0
        @test RenewalDiD._foi(x -> x == 2 ? 1 : 0, 1, 1, 10, 2) == 0.2
    end
    @testset "event rates" begin
        @test RenewalDiD._simulatedinfections(0, 50, 1, 1, 100, 1) == 0
        @test RenewalDiD._simulatedinfections(1, 50, 1, 1, 100, 1) == 1
        @test RenewalDiD._diseaseprogression(0, 0, 1) == 0
        @test RenewalDiD._diseaseprogression(10_000, 0.1, 1) == 1000
        @test RenewalDiD._diseaseprogression(10_000, x -> x == 1 ? 1 : 0, 1) == 10_000
        @test RenewalDiD._diseaseprogression(10_000, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._recovery(0, 0, 1) == 0
        @test RenewalDiD._recovery(500, 0.1, 1) == 50
        @test RenewalDiD._recovery(500, x -> x == 1 ? 1 : 0, 1) == 500
        @test RenewalDiD._recovery(500, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._recovery(100, 0.5, 1) == 50
        @test RenewalDiD._diagnosis(100, 0.5, 0, 1) == 0
        @test RenewalDiD._diagnosis(100, 0.5, 0.5, 1) == 50
        @test_throws ArgumentError RenewalDiD._diagnosis(100, 0.5, 1, 1)
        @test RenewalDiD._diagnosis(100, 1, 0.5, 1) == 100
        @test RenewalDiD._seirrates([0, 0, 0, 0, 1, 0], 1, 0, 0, 0, 0) == [0, 0, 0, 0, 0]
        @test_throws ArgumentError RenewalDiD._seirrates(
            [0, 0, 0, 0, 1, 0, 0], 1, 0, 0, 0, 0
        ) 
        @test RenewalDiD._seirrates(
            [10_000, 0, 500, 0, 9500, 0], 1, 0.5, 0, 0, 0
        ) == [125, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 0, 500, 500, 9000, 0], 1, 0.5, 0, 0, 0
        ) == [250, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 500, 500, 8000, 0], 1, 0.5, 0, 0, 0
        ) == [250, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 500, 500, 8000, 0], 1, 0.5, 0.5, 0, 0
        ) == [250, 500, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 1000, 0, 8000, 0], 1, 0.5, 0.5, 0.4, 0
        ) == [250, 500, 0, 400, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 1000, 1000, 7000, 0], 1, nothing, 0.5, 0.4, 0
        ) == [0, 500, 0, 400, 400]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 1000, 0, 8000, 0], 1, 0.5, 0.5, 0.4, 0.5
        ) == [250, 500, 400, 400, 0]
    end
    @testset "next event" begin
        v1 = 0.002143243061437003
        @test RenewalDiD._tstep(StableRNG(1), [100, 10, 50, 30, 60]) == v1
        @test RenewalDiD._nexteventtime(StableRNG(1), 3, [100, 10, 50, 30, 60]) == 3 + v1
        @test RenewalDiD._nexteventtime(3, [100, 10, 50, 30, 60]) > 3
        @test RenewalDiD._nextevent([100, 0, 0, 0, 0]) == 1
        @test RenewalDiD._nextevent(StableRNG(1), [0, 0, 100, 0, 0]) == 3
        @test RenewalDiD._nextevent(StableRNG(1), [100, 10, 50, 30, 60]) == 3
        @test RenewalDiD._nextevent(StableRNG(2), [100, 10, 50, 30, 60]) == 5
        u1 = [100, 100, 100, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 1)
        @test u1 == [99, 101, 100, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 2)
        @test u1 == [99, 100, 101, 100, 100, 0]
        RenewalDiD._updateevent!(u1, 3)
        @test u1 == [99, 100, 100, 101, 100, 1]
        RenewalDiD._updateevent!(u1, 4)
        @test u1 == [99, 100, 99, 101, 101, 1]
        RenewalDiD._updateevent!(u1, 5)
        @test u1 == [99, 100, 99, 100, 102, 1]
        u2 = [1, 0, 0, 1, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u2, 1, 0, 0, 0, 0)
        @test u2 == [1, 0, 0, 1, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u2, 1, 100, 0, 0, 0)
        @test u2 == [0, 1, 0, 1, 0, 0]
        u3 = [0, 1, 0, 0, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u3, 1, 0, 1000, 0, 0)
        @test u3 == [0, 0, 1, 0, 0, 0]
        u4 = [0, 0, 1, 0, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u4, 1, 0, 0, 1, 0.999)
        @test u4 == [0, 0, 0, 1, 0, 1]
        u5 = [0, 0, 1, 1, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u5, 1, 0, 0, 10, 0)
        @test u5 == [0, 0, 0, 0, 2, 0]
    end
    @testset "run whole simulation" begin
        @test_throws ArgumentError runsimulation(
            20, [100, 10, 50, 30, 60], 0.5, 0.5, 0.4, 0.5
        )
        @test_throws MethodError runsimulation(
            20.2, [100, 10, 50, 30, 60, 0], 0.5, 0.5, 0.4, 0.5
        )
        @test_throws ErrorException runsimulation(
            20, [100, 10, 50, 30, 60, 0], x -> 0.5 - 0.1 * x, 0.5, 0.4, 0.5
        )
        @test_throws ArgumentError runsimulation(
            20, [100, 10, 50, 30, 60, 0], 1.5, 0.5, 0.4, x -> 0.5 + 0.1 * x
        )
        S1 = runsimulation(StableRNG(1), 20, [100, 10, 50, 30, 60, 0], 0, 0.5, 0.4, 0.5)
        @test S1 == [
            100  10  50  30   60   0
            100   7  24  34   85  13
            100   3  14  25  108  17
            100   2   6  22  120  20
            100   1   4  16  129  22
            100   0   2   9  139  23
            100   0   0   7  143  24
            100   0   0   4  146  24
            100   0   0   2  148  24
            100   0   0   1  149  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
            100   0   0   0  150  24
        ]
        S2 = runsimulation(StableRNG(1), 20, [100, 10, 50, 30, 60, 0], 0.5, 0.5, 0.4, 0.5)
        @test S2 == [
            100  10  50  30   60   0
             90  16  23  33   88  16
             82  14  18  30  106  24
             72  17   8  23  130  31
             71  13   6  24  136  35
             68  11   5  16  150  37
             65  11   4  12  158  39
             62  11   5   6  166  40
             59   8   7   6  170  44
             56   7   8   3  176  45
             55   6   6   3  180  47
             55   4   3   4  184  50
             55   2   2   3  188  51
             55   1   2   3  189  51
             54   1   3   2  190  51
             54   1   2   0  193  52
             54   0   2   0  194  53
             53   1   0   1  195  54
             53   0   1   1  195  54
             53   0   1   0  196  54
             52   1   1   0  196  54
        ]
        @test simulationcases(zeros(10, 6)) == zeros(10)
        @test simulationcases(zeros(12, 6)) == zeros(12)
        @test_throws ArgumentError simulationcases(zeros(12, 4))
        SC1 = [
             0
            13
             4
             3
             2
             1
             1
             0
             0
             0
             0
             0
             0
             0
             0
             0
             0
             0
             0
             0
             0
        ]
        @test simulationcases(S1) == SC1
        @test simulationcases(
            StableRNG(1), 20, [100, 10, 50, 30, 60, 0], 0, 0.5, 0.4, 0.5
        ) == SC1
        SC2 = simulationcases(20, [100, 10, 50, 30, 60, 0], 0.5, 0.5, 0.4, 0.5)
        @test SC2 isa Vector{Int}
        @test length(SC2) == 21
    end 
end  # @testset "simulations"

@testset "fitting parameters" begin
    @test RenewalDiD._ntimes(zeros(2, 3)) == 2
    @test RenewalDiD._ntimes(zeros(3, 2)) == 3
    M1 = InterventionMatrix(4, [2, 3, nothing]) 
    @test RenewalDiD._ntimes(M1) == 4
    @test RenewalDiD._ngroups(zeros(2, 3)) == 3
    @test RenewalDiD._ngroups(zeros(3, 2)) == 2
    @test RenewalDiD._ngroups(M1) == 3
    @test_throws MethodError RenewalDiD._ntimes(zeros(2))
    @test_throws MethodError RenewalDiD._ngroups(zeros(2))
    @test RenewalDiD._gammavec(0, 1, zeros(3)) == zeros(3)
    @test RenewalDiD._gammavec(0, 1, zeros(4)) == zeros(4)
    @test RenewalDiD._gammavec(1, 1, zeros(3)) == ones(3)
    @test RenewalDiD._gammavec(1, 1, [1, -1, 0]) == [2, 0, 1]
    @test RenewalDiD._gammavec(1, 0.5, [1, -1, 0]) == [1.5, 0.5, 1]
    @test RenewalDiD._thetavec(0, zeros(3), 1) == zeros(4)
    @test RenewalDiD._thetavec(0, zeros(4), 1) == zeros(5)
    @test RenewalDiD._thetavec(1, zeros(3), 1) == ones(4)
    @test RenewalDiD._thetavec(1, ones(3), 1) == [1, 2, 3, 4]
    @test RenewalDiD._thetavec(1, ones(3), 0.5) == [1, 1.5, 2, 2.5]
    @test RenewalDiD._thetavec(1, [-2.5, 0.5, 0.5], 0.5) == [1, -0.25, 0, 0.25]
    @test RenewalDiD._predictedlogR_0(0, zeros(3), zeros(4), 0, M1) == zeros(4, 3)
    M2 = InterventionMatrix(5, [nothing, nothing]; mutewarnings=true) 
    @test RenewalDiD._predictedlogR_0(0, zeros(2), zeros(5), 2, M2) == zeros(5, 2)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(2), zeros(4), 0, M1)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(3), zeros(3), 0, M1)
    @test RenewalDiD._predictedlogR_0(1, zeros(3), zeros(4), 0, M1) == ones(4, 3)
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], zeros(4), 0, M1) == [
        0  1  2
        0  1  2
        0  1  2
        0  1  2
    ]
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], 0, M1) == [
        -1    0    1
         0.5  1.5  2.5
         2.5  3.5  4.5
         0    1    2
    ]
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], -1, M1) == [
        -1    0    1
        -0.5  1.5  2.5
         1.5  2.5  4.5
        -1    0    2
    ]
    @test RenewalDiD._expectedseedcases(zeros(20, 2), 7) == zeros(7, 2)
    @test RenewalDiD._expectedseedcases(zeros(20, 2), 5) == zeros(5, 2)
    @test RenewalDiD._expectedseedcases(zeros(20, 3), 5) == zeros(5, 3)
    obs1 = (_obs = zeros(Int, 20, 2); _obs[1, :] .+= 1; _obs[1:5, :] .+= 1; _obs)
    @test RenewalDiD._expectedseedcases(obs1, 5) == [
         0                          0
         0                          0
         0                          0
         0                          0
        (log(6 / 5) - log(2) / 5)  (log(6 / 5) - log(2) / 5)
    ]
    @test RenewalDiD._expectedseedcases(obs1, 5; doubletime=5) == [
         0                          0
         0                          0
         0                          0
         0                          0
        (log(6 / 5) - log(2) / 5)  (log(6 / 5) - log(2) / 5)
    ]
    @test RenewalDiD._expectedseedcases(obs1, 5; doubletime=10) == [
         0                               0
         0                               0
         0                               0
        (log(6 / 5) - 2 * log(2) / 10)  (log(6 / 5) - 2 * log(2) / 10) 
        (log(6 / 5) - log(2) / 10)      (log(6 / 5) - log(2) / 10)
    ]
    @test RenewalDiD._expectedseedcases(obs1, 5; sampletime=5) == [
         0                          0
         0                          0
         0                          0
         0                          0
        (log(6 / 5) - log(2) / 5)  (log(6 / 5) - log(2) / 5)
    ]
    @test RenewalDiD._expectedseedcases(obs1, 5; sampletime=3) == [
         0                              0
         0                              0
         0                              0
        (log(4 / 3) - 2 * log(2) / 5)  (log(4 / 3) - 2 * log(2) / 5)
        (log(4 / 3) - log(2) / 5)      (log(4 / 3) - log(2) / 5)
    ]
    obs2 = (_obs = zeros(Int, 20, 2); _obs[1, 1] += 1; _obs[1:5, :] .+= 1; _obs)
    @test RenewalDiD._expectedseedcases(obs2, 5) == [
         0                         0
         0                         0
         0                         0
         0                         0
        (log(6 / 5) - log(2) / 5)  0
    ]
    @test RenewalDiD._expectedseedcases(zeros(2, 3), 5) == zeros(5, 3)
    @test_throws BoundsError RenewalDiD._expectedseedcases(zeros(2, 3), 5; sampletime=5)
    # showing `log(0)` and `log(1)` to emphasize that these parameters are natural logarithms
    @test RenewalDiD._expectedinfections(g_covid, log(0), zeros(10)) == 0
    @test RenewalDiD._expectedinfections(g_covid, log(1), ones(10)) == 
        sum(RenewalDiD.COVIDSERIALINTERVAL[2:11])
    @test RenewalDiD._expectedinfections(g_covid, log(1), ones(6)) == 
        sum(RenewalDiD.COVIDSERIALINTERVAL[2:7])
    @test RenewalDiD._expectedinfections(g_covid, log(1), [0, 0, 0, 0, 0, 1]) == 0.0440204506
    @test RenewalDiD._expectedinfections(g_covid, log(1), [1, 0, 0, 0, 0, 0]) == 0.0917470443
    @test RenewalDiD._expectedinfections(g_seir, log(0), zeros(10); gamma=0.5) == 0
    @test RenewalDiD._expectedinfections(g_seir, log(1), [1, 0, 0, 0, 0, 0]; gamma=0.5) == 
        0.07468060255179591
    @test_throws UndefKeywordError RenewalDiD._expectedinfections(g_seir, log(0), zeros(10))
    # we do not need to use `generationtime` with `g_seir` but test it as though `g_seir` 
    # were a user-generated function
    @test RenewalDiD._expectedinfections(
        generationtime, log(0), zeros(10); 
        func=g_seir, gamma=0.5,
    ) == 0
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), [1, 0, 0, 0, 0, 0]; 
        func=g_seir, gamma=0.5,
    ) == 0.07468060255179591
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), [1, 0, 0, 0, 0, 0]; 
        func=g_seir, gamma=0.5, t_max=5,
    ) == 0
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), zeros(10); 
        vec=[0, 0, 0, 0, 0, 1],
    ) == 0
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), [1, 0, 0, 0, 0, 0]; 
        vec=[0, 0, 0, 0, 0, 0, 1],
    ) == 1
    @test RenewalDiD._expectedinfections(
        generationtime, log(2), ones(3); 
        vec=[0, 0, 0, 1, 0, 0, 0],
    ) == 2
    @test_throws ArgumentError RenewalDiD._expectedinfections(
        generationtime, log(0), zeros(10)
    ) 
    @test_throws ArgumentError RenewalDiD._expectedinfections(
        generationtime, log(0), zeros(10);
        func=g_seir, gamma=0.5, vec=[0, 0, 0, 0, 0, 1],
    ) 
    @test RenewalDiD._approxcases(0, 0) == 0
    @test RenewalDiD._approxcases(1, 0) == 1
    @test RenewalDiD._approxcases(1, 0.5) == 1.5
    @test RenewalDiD._approxcases(2, 0.5) == 3
    @test RenewalDiD._approxcases(2, -1.5) == 0
    @test RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(2), 2
    ) == zeros(5, 2)
    @test RenewalDiD._infections(
        g_covid, zeros(6, 2), log.(zeros(4, 2)), zeros(2, 2), zeros(2), 2
    ) == zeros(6, 2)
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 3)), zeros(2, 2), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 3), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(3), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(4, 2)), zeros(2, 2), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(2), 3
    ) 
    seedinfections1 = [
        1  0
        2  1
    ]
    @test RenewalDiD._infections(
        generationtime, zeros(5, 2), log.(zeros(3, 2)), seedinfections1, 1000 .* ones(2), 2;
        vec=[0, 1]
    ) == [
        1  0
        2  1 
        0  0
        0  0
        0  0
    ]
    @test RenewalDiD._infections(
        generationtime, 
        zeros(5, 2), 
        log.(1.5 .* ones(3, 2)), 
        seedinfections1, 
        1000 .* ones(2), 
        2;
        vec=[0, 1]
    ) == [
        1     0
        2     1 
        3     1.5
        4.5   2.25
        6.75  3.375
    ]
    M_x1 = [
         2    0
        -0.5  1
         1    0.5
        -2    1
         0    1
    ]
    @test RenewalDiD._infections(
        generationtime, M_x1, log.(ones(3, 2)), seedinfections1, 1000 .* ones(2), 2;
        vec=[0, 1]
    ) == [
        3   0
        1   2 
        2   3
        0   6
        0  12
    ]
end  # @testset "fitting parameters"

end  # @testset "RenewalDiD.jl"
