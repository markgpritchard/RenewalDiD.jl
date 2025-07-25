#package tests

# note that there are not yet tests of the `renewaldid` function itself

using RenewalDiD
using Test

@testset "RenewalDiD.jl" begin
    @testset "test functions for tests" begin
        include("functionsfortests_test.jl")
    end
    @testset "intervention matrix" begin
        include("interventionmatrix_test.jl")
    end 
    @testset "generation interval" begin
        include("generationinterval_test.jl")
    end 
    @testset "simulations" begin
        include("simulations_test.jl")
    end 
    @testset "fitting parameters (supplementary functions)" begin
        include("fittingparameters_supplementaryfunctionstest.jl")
    end  
    @testset "fitting parameters" begin
        include("fittingparameters_test.jl") 
    end    
    @testset "process parameters" begin
        include("processparameters_test.jl")
    end
    @testset "plotting without Makie" begin
        include("plotting_testwithoutMakie.jl")
    end
end 
