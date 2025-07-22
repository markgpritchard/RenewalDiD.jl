# test functions for processing output of fitted parameters

using RenewalDiD
using RenewalDiD.FittedParameterTestFunctions
using StableRNGs
using Test
using Turing

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
df2 = testdataframe( ; 
    nchains=4, 
    niterations=200, 
    ngroups=4, 
    ntimes=12, 
    nseeds=7,
)
df3 = testdataframe( ;
    nchains=1,
    niterations=1,
    ngroups=2,
    ntimes=2,
    nseeds=2,
    alpha=zeros(1),
    gammadefault=zeros(1), 
    thetadefault=zeros(1), 
    mxdefault=zeros(1),
)
df4 = testdataframe( ;
    nchains=1,
    niterations=1,
    ngroups=2,
    ntimes=2,
    nseeds=2,
    alpha=[log(2)],
    gammadefault=zeros(1), 
    thetadefault=zeros(1), 
    mxdefault=zeros(1),
)
df5 = let 
    kw = (Symbol("gammas_raw[1]") => [log(2)], )
    testdataframe( ;
        nchains=1,
        niterations=1,
        ngroups=2,
        ntimes=2,
        nseeds=2,
        alpha=zeros(1),
        sigma_gamma=ones(1),
        gammadefault=zeros(1), 
        thetadefault=zeros(1), 
        mxdefault=zeros(1),
        kw...
    )
end
df6 = let 
    kw = (Symbol("gammas_raw[1]") => [2 * log(2)], )
    testdataframe( ;
        nchains=1,
        niterations=1,
        ngroups=2,
        ntimes=2,
        nseeds=2,
        alpha=zeros(1),
        sigma_gamma=[0.5],
        gammadefault=zeros(1), 
        thetadefault=zeros(1), 
        mxdefault=zeros(1),
        kw...
    )
end
s1 = samplerenewaldidinfections(
    zeros(2), df1, zeros(10, 3), 100 .* ones(3), zeros(7, 3), 2, 3, 10
)
s2 = samplerenewaldidinfections(
    zeros(2), df2, zeros(12, 4), 100 .* ones(4), zeros(7, 4), 2, 4, 12
)
s3 = samplerenewaldidinfections(
    [0, 1], df3, zeros(2, 2), 100 .* ones(2), [0  0; 1  1], 1, 2, 2
) 
s4 = samplerenewaldidinfections(
    [0, 1], df4, zeros(2, 2), 100 .* ones(2), [0  0; 1  1], 1, 2, 2
)
s5 = samplerenewaldidinfections(
    [0, 1], df5, zeros(2, 2), 100 .* ones(2), [0  0; 1  1], 1, 2, 2
) 
s6 = samplerenewaldidinfections(
    [0, 1], df6, zeros(2, 2), 100 .* ones(2), [0  0; 1  1], 1, 2, 2
)
s7 = samplerenewaldidinfections(
    [0, 1], df6, 1;
    interventions=zeros(2, 2), 
    Ns=100 .* ones(2), 
    seedmatrix=[0  0; 1  1], 
    ngroups=2, 
    ntimes=2,
)

@testset "samples with no infections" begin
    @test s1 == zeros(11, 3)
    @test s2 == zeros(13, 4)
end

@testset "sampling errors" begin
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, zeros(10, 3), 100 .* ones(3), zeros(7, 4), 2, 3, 10
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, zeros(10, 3), 100 .* ones(3), zeros(8, 3), 2, 3, 10
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, zeros(10, 3), 100 .* ones(3), zeros(6, 3), 2, 3, 10
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, zeros(10, 2), 100 .* ones(2), zeros(7, 3), 2, 2, 10
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, zeros(10, 4), 100 .* ones(4), zeros(7, 3), 2, 4, 10
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, zeros(9, 3), 100 .* ones(3), zeros(7, 3), 2, 3, 9
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, zeros(11, 3), 100 .* ones(3), zeros(7, 3), 2, 3, 11
    ) 
    @test_throws BoundsError samplerenewaldidinfections(
        zeros(2), df1, zeros(10, 3), 100 .* ones(3), zeros(7, 3), 3, 3, 10
    )
end

@testset "samples with infections" begin
    @test s3 == [1  1; 1  1; 1  1]
    @test s4 == [1  1; 2  2; 4  4]
    @test s5 == [1  1; 2  0.5; 4  0.25]
    @test s6 == [1  1; 2  0.5; 4  0.25]
    @test s7 == [1  1; 2  0.5; 4  0.25]


end








