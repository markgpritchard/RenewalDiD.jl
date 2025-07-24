# test functions for processing output of fitted parameters

using RenewalDiD
using RenewalDiD.FittedParameterTestFunctions
using StableRNGs
using Test
using Turing

rng1 = StableRNG(1)

rankvaluedf8 = DataFrame(
    :iteration => repeat(1:1000; outer=4),
    :chain => repeat(1:4; inner=1000),
    :tau => [
        [rand() for _ in 1:1000]; 
        [1 + rand() for _ in 1:1000];
        [2 + rand() for _ in 1:1000];
        [3 + rand() for _ in 1:1000]
    ]
)
rankvaluedf9 = DataFrame(
    :iteration => repeat(1:1000; outer=4),
    :chain => repeat(1:4; inner=1000),
    :tau => [
        [rand() for _ in 1:1000]; 
        [2 + rand() for _ in 1:1000];
        [1 + rand() for _ in 1:1000];
        [3 + rand() for _ in 1:1000]
    ]
)
rankvaluedf10 = DataFrame(
    :iteration => repeat(1:100; outer=4),
    :chain => repeat(1:4; inner=100),
    :tau => [
        [rand() for _ in 1:100]; 
        [1 + rand() for _ in 1:100];
        [2 + rand() for _ in 1:100];
        [3 + rand() for _ in 1:100]
    ]
)
rankvaluedf11 = DataFrame(
    :iteration => repeat(1:100; outer=3),
    :chain => repeat(1:3; inner=100),
    :tau => [
        [rand() for _ in 1:100]; 
        [1 + rand() for _ in 1:100];
        [2 + rand() for _ in 1:100]
    ]
)
rankvaluedf12 = DataFrame(
    :iteration => repeat(1:100; outer=3),
    :chain => repeat([1, 2, 4]; inner=100),
    :tau => [
        [rand() for _ in 1:100]; 
        [1 + rand() for _ in 1:100];
        [2 + rand() for _ in 1:100]
    ]
)
rankvaluedf13 = DataFrame(
    :iteration => repeat(1:1000; outer=4),
    :chain => repeat([1, 2, 2, 3]; inner=1000),
    :tau => rand(4000)
)
rankvaluedf14 = DataFrame(
    :iteration => repeat(51:150; outer=3),
    :chain => repeat([1, 2, 4]; inner=100),
    :tau => [
        [rand() for _ in 1:100]; 
        [1 + rand() for _ in 1:100];
        [2 + rand() for _ in 1:100]
    ]
)
rankvaluedf15 = DataFrame(
    :iteration => [1:9; 1:10; 1:10],
    :chain => [[1 for _ in 1:9]; [2 for _ in 1:10]; [3 for _ in 1:10]],
    :tau => rand(29)
)
rankvaluedf16 = DataFrame(
    :iteration => repeat(1:10; outer=3),
    :chain => repeat(1:3; inner=10),
    :tau => [
        [rand() for _ in 1:5]; 
        [1 + rand() for _ in 1:5];
        [1 + rand() for _ in 1:5]; 
        [2 + rand() for _ in 1:5];
        [2 + rand() for _ in 1:5]; 
        [rand() for _ in 1:5];
    ]
)
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
df7 = let 
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
        mxdefault=(-1.5 .* ones(1)),
        kw...
    )
end

s1 = samplerenewaldidinfections(
    zeros(2), df1, 2;
    interventions=zeros(10, 3), 
    Ns=(100 .* ones(3)), 
    seedmatrix=zeros(7, 3), 
    ngroups=3, 
    ntimes=10,
)
s2 = samplerenewaldidinfections(
    zeros(2), df2, 2;
    interventions=zeros(12, 4), 
    Ns=(100 .* ones(4)), 
    seedmatrix=zeros(7, 4), 
    ngroups=4, 
    ntimes=12,
)
s3 = samplerenewaldidinfections(
    [0, 1], df3, 1;
    interventions=zeros(2, 2), 
    Ns=(100 .* ones(2)), 
    seedmatrix=[0  0; 1  1], 
    ngroups=2, 
    ntimes=2,
)
s4 = samplerenewaldidinfections(
    [0, 1], df4, 1;
    interventions=zeros(2, 2), 
    Ns=(100 .* ones(2)), 
    seedmatrix=[0  0; 1  1], 
    ngroups=2, 
    ntimes=2,
)
s5 = samplerenewaldidinfections(
    [0, 1], df5, 1;
    interventions=zeros(2, 2), 
    Ns=(100 .* ones(2)), 
    seedmatrix=[0  0; 1  1], 
    ngroups=2, 
    ntimes=2,
)
s6 = samplerenewaldidinfections(
    [0, 1], df6, 1;
    interventions=zeros(2, 2), 
    Ns=(100 .* ones(2)), 
    seedmatrix=[0  0; 1  1], 
    ngroups=2, 
    ntimes=2,
)
s7 = samplerenewaldidinfections(
    [0, 1], df7, 1;
    interventions=zeros(2, 2), 
    Ns=(100 .* ones(2)), 
    seedmatrix=[0  0; 1  1], 
    ngroups=2, 
    ntimes=2,
)

rv16a = let
    _rvs = rankvalues(rankvaluedf16, :tau)
    _a = _rvs[1]
    _b = _rvs[11]
    _c = _rvs[21]
    repeat([_a, _b, _c]; inner=10)
end

rv16b = let
    _rvs = rankvalues(rankvaluedf16, :tau; binsize=7)
    _a1 = _rvs[1]
    _a2 = _rvs[8]
    _b1 = _rvs[11]
    _b2 = _rvs[18]
    _c1 = _rvs[21]
    _c2 = _rvs[28]
    [
        [_a1 for _ in 1:7]; 
        [_a2 for _ in 1:3]; 
        [_b1 for _ in 1:7]; 
        [_b2 for _ in 1:3]; 
        [_c1 for _ in 1:7]; 
        [_c2 for _ in 1:3]; 
    ]
end

sim1 = testsimulation(rng1)
model1 =  renewaldid(
    sim1, g_seir, packpriors(; sigma_thetaprior=Exponential(0.05)); 
    gamma=0.2, sigma=0.5
)
chain1 = sample(rng1, model1, NUTS(), MCMCThreads(), 20, 4; verbose=false, progress=false)
df7 = DataFrame(chain1)

@testset "number of unique elements" begin
    @test nunique([1, 2, 3]) == 3
    @test nunique(ones(3)) == 1
end

@testset "rank values for trace plots" begin
    @test rankvalues(rankvaluedf8, :tau) == repeat(1:4; inner=1000) 
    @test rankvalues(rankvaluedf9, :tau) == repeat([1, 3, 2, 4]; inner=1000) 
    @test rankvalues(rankvaluedf10, :tau) == repeat(1:4; inner=100) 
    @test rankvalues(rankvaluedf11, :tau) == repeat(1:3; inner=100) 
    @test rankvalues(rankvaluedf12, :tau) == repeat(1:3; inner=100)
    @test_throws ErrorException rankvalues(rankvaluedf13, :tau)
    @test rankvalues(rankvaluedf14, :tau) == repeat(1:3; inner=100)
    @test_throws ErrorException rankvalues(rankvaluedf15, :tau)
    @test rankvalues(rankvaluedf16, :tau; binsize=5) == repeat([1, 2, 2, 3, 3, 1]; inner=5)
    @test rankvalues(rankvaluedf16, :tau) == rv16a
    @test rankvalues(rankvaluedf16, :tau; binsize=7) == rv16b
end

@testset "samples with no infections" begin
    @test s1 == zeros(11, 3)
    @test s2 == zeros(13, 4)
end

@testset "samples with infections, M_x all 0" begin
    @test s3 == [1  1; 1  1; 1  1]
    @test s4 == [1  1; 2  2; 4  4]
    @test s5 == [1  1; 2  0.5; 4  0.25]
    @test s6 == [1  1; 2  0.5; 4  0.25]
end

@testset "samples with Mx < -1" begin
    @test s7 == zeros(3, 2)

end

@testset "keyword errors" begin
    @test_throws ArgumentError samplerenewaldidinfections(zeros(2), df1, 2)
    @test_throws ArgumentError samplerenewaldidinfections(
        g_seir, df7, 1; 
        data=sim1, gamma=0.2, sigma=0.5,
    )
end

@testset "sampling errors" begin
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, 2;
        interventions=zeros(10, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 4),
        ngroups=3,
        ntimes=10,
    )
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, 2;
        interventions=zeros(10, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(8, 3), 
        ngroups=3, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, 2;
        interventions=zeros(10, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(6, 3), 
        ngroups=3, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, 2;
        interventions=zeros(10, 2), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 3), 
        ngroups=2, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, 2;
        interventions=zeros(10, 4), 
        Ns=100 .* ones(4), 
        seedmatrix=zeros(7, 3), 
        ngroups=4, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, 2;
        interventions=zeros(9, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 3), 
        ngroups=3, 
        ntimes=9,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1, 2;
        interventions=zeros(11, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 3), 
        ngroups=3, 
        ntimes=11,
    ) 
    @test_throws BoundsError samplerenewaldidinfections(
        zeros(2), df1, 3;
        interventions=zeros(10, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 3), 
        ngroups=3, 
        ntimes=10,
    ) 
end
