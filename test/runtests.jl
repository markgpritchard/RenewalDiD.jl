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
    @test isnothing(testgenerationtime([0.0, 0.2, 0.4]; muteinfo=true))
    @test generationtime([0.0, 0.2, 0.4], 0) == 0
    @test generationtime([0.0, 0.2, 0.4], 1) == 0.2
    @test generationtime([0.0, 0.2, 0.4], -1) == 0
    @test generationtime([0.0, 0.2, 0.4], 3) == 0
    @test_throws MethodError generationtime([0.0, 0.2, 0.4], 2.3)
    #@test_throws ArgumentError generationtime([0.0, 0.2, 0.4, 0.8], 0)
    # this function no longer throws an error in this case. Use:
    @test_throws ArgumentError testgenerationtime([0.0, 0.2, 0.4, 0.8])
    #@test_throws ArgumentError generationtime([0.0, 0.2, 0.4, -0.1], 0)
    # this function no longer throws an error in this case. Use:
    @test_throws ArgumentError testgenerationtime([0.0, 0.2, 0.4, -0.1])
    @test generationtime(x -> x < 2 ? 0.0 : 1/x^2, 0) == 0
    @test generationtime(x -> x < 2 ? 0.0 : 1/x^2, 2) == 0.25
    @test generationtime(x -> x == -1 ? 0.25 : 0, -1) == 0
    @test_throws MethodError generationtime(x -> x < 2 ? 0.0 : 1/x^2, 1.5)
    #@test_throws ArgumentError generationtime(x -> 1/x^2, 0)
    # this function no longer throws an error in this case. Use:
    @test_throws ArgumentError testgenerationtime(x -> 1/x^2)
    #@test_throws ArgumentError generationtime(x -> x == 12 ? -0.5 : 0.0, 0)
    # this function no longer throws an error in this case. Use:
    @test_throws ArgumentError testgenerationtime(x -> x == 12 ? -0.5 : 0.0)
    @test generationtime(x -> x == 12 ? -0.5 : 0.0, 1; t_max=10) == 0
    @test generationtime(x -> x == 12 ? -0.5 : 0.0, 12; t_max=10) == 0
    @test isnothing(testgenerationtime(x -> x == 12 ? -0.5 : 0.0; t_max=10, muteinfo=true))
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
    @test isnothing(testgenerationtime(g_covid; muteinfo=true))
    @test isnothing(testgenerationtime(g_seir; gamma=0.4, sigma=0.5, muteinfo=true))
    @test isnothing(testgenerationtime(g_seir; gamma=0.4, muteinfo=true))
    @test isnothing(testgenerationtime(vectorg_seir(0.4); muteinfo=true))
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
        @test RenewalDiD._parameter(x -> 2 * x, 0; upper=1) == 0 
        @test RenewalDiD._parameter(x -> 2 * x, 0.5; upper=1) == 1 
        @test_throws ErrorException RenewalDiD._parameter(x -> 2 * x, 1; upper=1)
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
        @test RenewalDiD._simulatedinfections(0, 50, 1, 1, 1, 100, 1) == 0
        @test RenewalDiD._simulatedinfections(1, 50, 1, 1, 1, 100, 1) == 1.5
        @test RenewalDiD._diseaseprogression_to_i_n(0, 0, 0, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0, 0.1, 1) == 1000
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 1, 0.1, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0.25, 0.1, 1) == 750
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_n(1000, 1.5, 0.1, 1)
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, x -> x == 1 ? 1 : 0, 1, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, x -> x == 1 ? 1 : 0, 1, 2) == 
            10_000
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, x -> x == 2 ? 1.5 : 0, 1, 1) == 
            10_000
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_n(
            10_000, x -> x == 2 ? 1.5 : 0, 1, 2
        )
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0, x -> x == 1 ? 1 : 0, 1) == 
            10_000
        @test RenewalDiD._diseaseprogression_to_i_n(10_000, 0, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._diseaseprogression_to_i_f(0, 0, 0, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 0, 0.1, 1) == 0
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 1, 0.1, 1) == 1000
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 0.25, 0.1, 1) == 250
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_f(1000, 1.5, 0.1, 1)
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, x -> x == 1 ? 1 : 0, 1, 1) == 
            10_000
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, x -> x == 1 ? 1 : 0, 1, 2) == 0
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, x -> x == 2 ? 1.5 : 0, 1, 1) == 0
        @test_throws ErrorException RenewalDiD._diseaseprogression_to_i_f(
            10_000, x -> x == 2 ? 1.5 : 0, 1, 2
        )
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 1, x -> x == 1 ? 1 : 0, 1) == 
            10_000
        @test RenewalDiD._diseaseprogression_to_i_f(10_000, 1, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._diagnosis(0, 0, 1) == 0
        @test RenewalDiD._diagnosis(500, 0.1, 1) == 50
        @test RenewalDiD._diagnosis(500, x -> x == 1 ? 1 : 0, 1) == 500
        @test RenewalDiD._diagnosis(500, x -> x == 1 ? 1 : 0, 2) == 0
        @test RenewalDiD._diagnosis(100, 0.5, 1) == 50
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
        @test RenewalDiD._n_seir(zeros(Int, 7)) == 0
        @test RenewalDiD._n_seir(ones(Int, 7)) == 6
        @test_throws ArgumentError RenewalDiD._n_seir(ones(Int, 8))
        @test_throws MethodError RenewalDiD._n_seir(ones(7))
        @test_throws ArgumentError RenewalDiD._n_seir([1, 1, 0, 0, 0, -1, 0])
        @test RenewalDiD._seirrates([0, 0, 0, 0, 0, 1, 0], 1, 0, 0, 0, 0, 0) == zeros(6)
        @test_throws ArgumentError RenewalDiD._seirrates(
            [0, 0, 0, 0, 1, 0], 1, 0, 0, 0, 0, 0
        ) 
        @test RenewalDiD._seirrates([10_000, 0, 500, 0, 0, 9500, 0], 1, 0.5, 0, 0, 0, 0) == 
            [125, 0, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 0, 500, 500, 0, 9000, 0], 1, 0.5, 0, 0, 0, 0
        ) == [250, 0, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 0, 500, 500, 500, 8500, 0], 1, 0.5, 0, 0, 0, 0
        ) == [375, 0, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 500, 500, 0, 8000, 0], 1, 0.5, 0, 0, 0, 0
        ) == [250, 0, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 500, 500, 0, 8000, 0], 1, 0.5, 0, 0, 0, 0.5
        ) == [250, 500, 0, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 500, 500, 0, 8000, 0], 1, 0.5, 0, 0, 1, 0.5
        ) == [250, 0, 500, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 1000, 0, 0, 8000, 0], 1, 0.5, 0, 0.5, 0.6, 0.5
        ) == [250, 200, 300, 0, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 600, 400, 0, 8000, 0], 1, 0.5, 0, 0.5, 0.6, 0.5
        ) == [250, 200, 300, 200, 0, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 600, 400, 0, 8000, 0], 1, 0.5, 1/3, 0.5, 0.6, 0.5
        ) == [250, 200, 300, 200, 200, 0]
        @test RenewalDiD._seirrates(
            [10_000, 1000, 600, 400, 1000, 7000, 0], 1, 0.5, 1/3, 0.5, 0.6, 0.5
        ) == [500, 200, 300, 200, 200, 1000]
    end
    @testset "next event" begin
        v1 = 0.002143243061437003
        @test RenewalDiD._tstep(StableRNG(1), [80, 10, 50, 30, 60, 20]) == v1
        @test RenewalDiD._nexteventtime(StableRNG(1), 3, [80, 10, 50, 30, 60, 20]) == 3 + v1
        @test RenewalDiD._nexteventtime(3, [80, 10, 50, 30, 60, 20]) > 3
        @test RenewalDiD._nextevent([100, 0, 0, 0, 0, 0]) == 1
        @test RenewalDiD._nextevent(StableRNG(1), [0, 0, 100, 0, 0, 0]) == 3
        @test RenewalDiD._nextevent(StableRNG(1), [100, 10, 50, 30, 60, 0]) == 3
        @test RenewalDiD._nextevent(StableRNG(2), [100, 10, 50, 30, 60, 0]) == 5
        @testset "update event" begin
            # this test set mutates the same vector multiple times so if multiple tests fail
            # it may only be a problem with the first of the failing functions
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
        # some functions in this test set mutate the same vector multiple times so if 
        # multiple tests fail it may only be a problem with the first of the failing 
        # functions
        u2 = [1, 0, 0, 1, 0, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u2, 1, 0, 0, 0, 0, 0)
        @test u2 == [1, 0, 0, 1, 0, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u2, 1, 100, 0, 0, 0, 0)
        @test u2 == [0, 1, 0, 1, 0, 0, 0]
        u3 = [0, 1, 0, 0, 0, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u3, 1, 0, 0, 0, 0, 100)
        @test u3 == [0, 0, 1, 0, 0, 0, 0]
        RenewalDiD._simulateday!(StableRNG(1), u3, 1, 0, 100, 1000, 0, 0)
        @test u3 == [0, 0, 0, 0, 0, 1, 0]
        u4 = [0, 1, 0, 0, 0, 0, 0]
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
        @test_throws MethodError runsimulation(
            20.2, [100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.4, 0.5
        )
        @test_throws ErrorException runsimulation(
            20, [100, 10, 50, 30, 40, 20, 0], x -> 0.5 - 0.1 * x, 0.4, 0.5, 0.4, 0.5
        )
        @test_throws ErrorException runsimulation(
            20, [100, 10, 50, 30, 40, 20, 0], 1.5, 0.4, 0.5, x -> 0.5 + 0.1 * x, 0.5
        )
        @test_throws ArgumentError runsimulation(
            -20, [100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5
        )
        S1 = runsimulation(
            StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5
        )
        @test S1 == [
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
        S2 = runsimulation(
            StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5
        )
        @test S2 == [
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
        @test simulationcases(zeros(10, 7)) == zeros(10)
        @test simulationcases(zeros(12, 7)) == zeros(12)
        @test_throws ArgumentError simulationcases(zeros(12, 6))
        SC1 = [
             0
            11
             7
             6
             2
             2
             0
             2
             1
             0
             1
             0
             2
             0
             0
             0
             0
             1
             0
             0
             0
        ]
        @test simulationcases(S1) == SC1
        @test simulationcases(
            StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5
        ) == SC1
        SC2 = simulationcases(
            StableRNG(1), 20, [100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5
        )
        @test SC2 isa Vector{Int}
        @test length(SC2) == 21
    end 
    @testset "group simulations" begin
        rng1 = StableRNG(1)
        sim1 = simulationcases(
            rng1, 20, [100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5
        )
        sim2 = simulationcases(
            rng1, 20, [100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5
        )
        sim3 = simulationcases(
            rng1, 
            20, 
            [100, 10, 50, 1, 40, 20, 0], 
            x -> max(0.1, 2 - 0.1 * x), 
            0.4, 
            0.5, 
            0.5, 
            0.5
        )
        rng2 = StableRNG(1)
        dict1 = packsimulations(
            rng2,
            20,
            ([100, 10, 50, 30, 40, 20, 0], 0, 0.4, 0.5, 0.5, 0.5, 10),
            ([100, 10, 50, 30, 40, 20, 0], 0.5, 0.4, 0.5, 0.5, 0.5, 12),
            ([100, 10, 50, 1, 40, 20, 0], x -> max(0.1, 2 - 0.1 * x), 0.4, 0.5, 0.5, 0.5, nothing),
        )
        @test dict1[:interventions] == InterventionMatrix(20, 10, 12, nothing)
        @test dict1[:Ns] == [250, 250, 221]
        @testset for i in 1:3 
            @test dict1[:observedcases][:, i] == [sim1, sim2, sim3][i]
        end
    end
end  # @testset "simulations"

@testset "fitting parameters" begin
    @testset "functions called by `RenewalDiD._renewaldid`" begin
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
        @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(
            0, zeros(2), zeros(4), 0, M1
        )
        @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(
            0, zeros(3), zeros(3), 0, M1
        )
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
        @test RenewalDiD._expectedinfections(g_covid, log(1), [0, 0, 0, 0, 0, 1]) ==    
            0.0440204506
        @test RenewalDiD._expectedinfections(g_covid, log(1), [1, 0, 0, 0, 0, 0]) == 
            0.0917470443
        @test RenewalDiD._expectedinfections(g_seir, log(0), zeros(10); gamma=0.5) == 0
        @test RenewalDiD._expectedinfections(
            g_seir, log(1), [1, 0, 0, 0, 0, 0]; 
            gamma=0.5
        ) == 0.07468060255179591
        @test_throws UndefKeywordError RenewalDiD._expectedinfections(
            g_seir, log(0), zeros(10)
        )
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
            generationtime, 
            zeros(5, 2), 
            log.(zeros(3, 2)), 
            seedinfections1, 
            1000 .* ones(2), 
            2;
            vec=[0, 1],
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
            vec=[0, 1],
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
    end
end  # @testset "fitting parameters"

end  # @testset "RenewalDiD.jl"
