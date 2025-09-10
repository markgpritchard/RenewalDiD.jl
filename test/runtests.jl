#package tests

using RenewalDiD
using Test

@testset "RenewalDiD.jl" begin
    @testset "documentation tests" begin
        #@info "test docstrings"
        #include("documenterdoctest.jl")
        @warn "not currently testing documenterdoctest.jl"
    end
    @testset "Aqua tests" begin
        @info "Aqua tests"
        include("aquatests.jl")
    end
    #=
    # there are no tests in functionsfortests_test.jl
    @testset "test functions for tests" begin
        include("functionsfortests_test.jl")
    end
    =#
    @testset "intervention vector" begin
        @info "starting intervention vector tests"
        include("interventionvector_test.jl")
    end 
    @testset "intervention matrix" begin
        @info "starting intervention matrix tests"
        include("interventionmatrix_test.jl")
    end 
    @testset "intervention array" begin
        @info "starting intervention array-3 tests"
        include("interventionarray3_test.jl")
    end 
    @testset "generation interval" begin
        @info "starting generation interval tests"
        include("generationinterval_test.jl")
    end 
    @testset "simulations" begin
        @info "starting simulations tests"
        include("simulations_test.jl")
    end 
    @testset "fitting parameters (supplementary functions)" begin
        @info "starting fitting parameters (supplementary functions) tests"
        include("fittingparameters_supplementaryfunctionstest.jl")
    end  
    @testset "fitting parameters" begin
        @info "starting fitting parameters"
        include("fittingparameters_test.jl") 
    end    
    @testset "process parameters" begin
        @info "starting process parameters tests"
        include("processparameters_test.jl")
    end
    @testset "plotting without Makie" begin
        @info "starting plotting without Makie tests"
        include("plotting_testwithoutMakie.jl")
    end
    @testset "example workflow" begin
        #@info "starting example workflow"
        #include("exampleworkflow.jl")
        @warn "not currently testing exampleworkflow.jl"
    end
end 
