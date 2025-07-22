# test generation interval functions

using RenewalDiD
using Test

v1 = [g_covid(x) for x in -1000:1:1000]
v2 = [g_seir(x; gamma=0.5) for x in -1000:1:1000]
v3 = [g_seir(x; gamma=0.4, sigma=0.5) for x in -1000:1:1000]
v4 = vectorg_seir(0.4, 0.5)
v5 = vectorg_seir(0.4, 0.5; t_max=10)
v6 = vectorg_seir(0.4)
v7 = vectorg_seir(0.4; t_max=10)

tw1 = "g(0) == 0.1: g(0) is never called and is assumed to equal 0. If you intended this \
    value for g(1) you should recheck the indexing"
tw2 = "g(0) == 0.2: g(0) is never called and is assumed to equal 0. If you intended this \
    value for g(1) you should recheck the indexing"

@testset "function `g_covid`" begin
    @testset "`g(0) == 0" begin
        @test g_covid(0) == 0
    end
    @testset "check specific values" begin
        @test g_covid(1) == 0.0440204506
        @test g_covid(30) == 0.0002057625
    end
    @testset "values outside range of the vector" begin
        @test g_covid(31) == 0
        @test g_covid(-1) == 0
    end
    @testset "`t` must be integer" begin
        @test_throws MethodError g_covid(1.2)
    end
    @testset "no negative values" begin
        @test minimum(v1) >= 0
    end
    @testset "sum of all values ≤ 1" begin
        @test sum(v1) <= 1
    end
    @testset "passes `testgenerationtime`" begin
        @test isnothing(testgenerationtime(g_covid; muteinfo=true))
    end
end

@testset "function `g_seir`" begin
    @testset "`g(0) == 0" begin
        @test g_seir(0; gamma=0.5) == 0
        @test g_seir(0; gamma=rand()) == 0
        @test g_seir(0; gamma=0.4, sigma=0.5) == 0
        @test g_seir(0; gamma=rand(), sigma=rand()) == 0
    end
    @testset "check specific values" begin
        @test g_seir(1; gamma=0.4, sigma=0.5) == 0.12757877264601186
        @test g_seir(1; gamma=0.5, sigma=0.5) == 0.15163266492815836
        @test g_seir(1; gamma=0.5) == 0.15163266492815836
        
    end
    @testset "cannot work without `gamma` keyword" begin
        @test_throws UndefKeywordError g_seir(1)
    end
    @testset "negative times" begin
        @test g_seir(-1; gamma=0.5) == 0
        @test g_seir(-1; gamma=0.4, sigma=0.5) == 0
    end
    @testset "no negative values" begin
        @test minimum(v2) >= 0
        @test minimum(v3) >= 0
    end
    @testset "sum of all values ≤ 1" begin
        @test sum(v2) <= 1
        @test sum(v3) <= 1
    end
    @testset "passes `testgenerationtime`" begin
        @test isnothing(testgenerationtime(g_seir; gamma=0.4, sigma=0.5, muteinfo=true))
        @test isnothing(testgenerationtime(g_seir; gamma=0.4, muteinfo=true))
    end
end

@testset "function `vectorg_seir`" begin
    @testset "equivalence to `g_seir`" begin
        @test generationtime(v4, 1) == g_seir(1; gamma=0.4, sigma=0.5)
        @test generationtime(v4, 12) == g_seir(12; gamma=0.4, sigma=0.5)
        @test generationtime(v4, 30) == 0
        @test generationtime(v5, 2) == g_seir(2; gamma=0.4, sigma=0.5)
        @test generationtime(v5, 12) == 0
        @test generationtime(v6, 1) == g_seir(1; gamma=0.4)
        @test generationtime(v6, 12) == g_seir(12; gamma=0.4, sigma=0.4)
        @test generationtime(v6, 30) == 0
        @test generationtime(v7, 2) == g_seir(2; gamma=0.4)
        @test generationtime(v7, 12) == 0
    end
    @testset "passes `testgenerationtime`" begin
        @test isnothing(testgenerationtime(vectorg_seir(0.4); muteinfo=true))
    end
end

@testset "user-specified vector generation intervals" begin
    @testset "check specific values" begin
        @test generationtime([0.0, 0.2, 0.4], 0) == 0
        @test generationtime([0.0, 0.2, 0.4], 1) == 0.2     
    end
    @testset "values outside range of the vector" begin
        @test generationtime([0.0, 0.2, 0.4], -1) == 0
        @test generationtime([0.0, 0.2, 0.4], 3) == 0
    end
    @testset "`t` must be integer" begin
        @test_throws MethodError generationtime([0.0, 0.2, 0.4], 2.3)
    end
    @testset "test generation time vectors with no issues" begin
        @test isnothing(testgenerationtime([0.0, 0.2, 0.4]; muteinfo=true))
    end
    @testset "test generation time with negative values" begin
        @test_throws ArgumentError testgenerationtime([0.0, 0.2, 0.4, -0.1])
    end
    @testset "test generation time that sums to >1" begin
        @test_throws ArgumentError testgenerationtime([0.0, 0.2, 0.4, 0.8])
    end
    @testset "warnings for g(0) != 0" begin
        @test_warn tw1 testgenerationtime([0.1, 0.2, 0.4]; muteinfo=true)
        @test_warn tw2 testgenerationtime([0.2, 0.2, 0.4]; muteinfo=true)
    end
    @testset "no warnings for g(0) == 0" begin
        @test_nowarn testgenerationtime([0.0, 0.2, 0.4]; muteinfo=true)
    end
end

@testset "user-specified function generation intervals" begin
    @testset "check specific values" begin
        @test generationtime(x -> x < 2 ? 0.0 : 1/x^2, 0) == 0
        @test generationtime(x -> x < 2 ? 0.0 : 1/x^2, 2) == 0.25
    end
    @testset "negative times" begin
        @test generationtime(x -> x == -1 ? 0.25 : 0, -1) == 0
    end
    @testset "`t` must be integer" begin
        @test_throws MethodError generationtime(x -> x < 2 ? 0.0 : 1/x^2, 1.5)
    end
    @testset "test generation time vectors with no issues" begin
        @test isnothing(testgenerationtime(x -> x == 5 ? 0.1 : 0; muteinfo=true))
    end
    @testset "test generation time with negative values" begin
        @test_throws ArgumentError testgenerationtime(x -> x == 12 ? -0.5 : 0.0)
    end
    @testset "test generation time that sums to >1" begin
        @test_throws ArgumentError testgenerationtime(x -> 1/x^2)
    end
    @testset "test generation time vectors with issues removed by `t_max`" begin
        @test generationtime(x -> x == 12 ? -0.5 : 0.0, 1; t_max=10) == 0
        @test generationtime(x -> x == 12 ? -0.5 : 0.0, 12; t_max=10) == 0
        @test isnothing(
            testgenerationtime(x -> x == 12 ? -0.5 : 0.0; t_max=10, muteinfo=true)
        )
    end
        @testset "warnings for g(0) != 0" begin
        @test_warn tw1 testgenerationtime(x -> x == 0 ? 0.1 : 0; muteinfo=true)
        @test_warn tw2 testgenerationtime(x -> x == 0 ? 0.2 : 0; muteinfo=true)
    end
    @testset "no warnings for g(0) == 0" begin
        @test_nowarn testgenerationtime(x -> x == 5 ? 0.1 : 0; muteinfo=true)
    end
end

@testset "correct functions for each generation function" begin
    let 
        _f(::RenewalDiD._Useablegenerationfunctions) = 1
        _f(x) = 2
        @test _f(zeros(Int, 2)) == 2
        @test _f(x -> 2 * x) == 2 
        @test _f(g_covid) == 1 
        @test _f(g_seir) == 1 
        @test _f(generationtime) == 1
    end    
end
