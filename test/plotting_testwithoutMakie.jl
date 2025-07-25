# tests of Plotting module without loading Makie 

using DataFrames
using RenewalDiD
using RenewalDiD.FittedParameterTestFunctions
using RenewalDiD.Plotting
using Test

df1 = testdataframe( ; 
    nchains=1, 
    niterations=2, 
    ngroups=3, 
    ntimes=10, 
    nseeds=7, 
    tau=[0, rand()], 
    alpha=[0, rand()], 
    sigma_gamma=[0, rand()], 
    sigma_theta=[0, rand()], 
)

@testset "method errors" begin
    @test_throws MethodError traceplot(df1, :tau)
    @test_throws MethodError tracerankplot(df1, :tau)
    @test_throws MethodError plotmodeloutput(df1, :tau)
end
