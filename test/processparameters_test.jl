# test functions for processing output of fitted parameters

using RenewalDiD
using RenewalDiD.FittedParameterTestFunctions
using StableRNGs
using Test
using Turing

rng1 = StableRNG(1)
_r1 = rand()
_r2 = rand()

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
    niterations=10,
    ngroups=2,
    ntimes=2,
    nseeds=2,
    alpha=zeros(10),
    gammadefault=zeros(10), 
    thetadefault=zeros(10), 
    mxdefault=zeros(10),
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
df8 = testdataframe( ;
    nchains=2,
    niterations=8,
    ngroups=3,
    ntimes=10,
    nseeds=2,
    alpha=zeros(16),
    gammadefault=zeros(16), 
    thetadefault=zeros(16), 
    mxdefault=zeros(16),
)
df9 = testdataframe( ;
    nchains=3,
    niterations=8,
    ngroups=3,
    ntimes=10,
    nseeds=2,
    alpha=zeros(24),
    gammadefault=zeros(24), 
    thetadefault=zeros(24), 
    mxdefault=zeros(24),
)
df10 = testdataframe( ;
    nchains=2,
    niterations=8,
    ngroups=4,
    ntimes=20,
    nseeds=2,
    alpha=zeros(16),
    gammadefault=zeros(16), 
    thetadefault=zeros(16), 
    mxdefault=zeros(16),
)
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
s8 = samplerenewaldidinfections(
    zeros(2), df8;
    interventions=zeros(10, 3), 
    Ns=(100 .* ones(3)), 
    seedmatrix=zeros(2, 3), 
    ngroups=3, 
    ntimes=10,
)
s9 = samplerenewaldidinfections(
    zeros(2), df9;
    interventions=zeros(10, 3), 
    Ns=(100 .* ones(3)), 
    seedmatrix=zeros(2, 3), 
    ngroups=3, 
    ntimes=10,
)
s10 = samplerenewaldidinfections(
    zeros(2), df10;
    interventions=zeros(20, 4), 
    Ns=(100 .* ones(4)), 
    seedmatrix=zeros(2, 4), 
    ngroups=4, 
    ntimes=20,
)
s3ma = samplerenewaldidinfections(
    zeros(2), df3;
    interventions=zeros(2, 2), 
    Ns=(100 .* ones(2)), 
    seedmatrix=[0  0; _r1  _r2], 
    ngroups=2, 
    ntimes=2,
)
s3mb = samplerenewaldidinfections(
    [0, 1], df3;
    interventions=zeros(2, 2), 
    Ns=(100 .* ones(2)), 
    seedmatrix=[0  0; 1  1], 
    ngroups=2, 
    ntimes=2,
)
s3mc = samplerenewaldidinfections(
    [0, 1], df3, 4:6;
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
chaindf1 = DataFrame(chain1)

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
    # test removed as not providing `n_seeds` to `samplerenewaldidinfections` no longer 
    # intended to throw an error
end


@testset "sample multiple rows" begin
    @test s8 == zeros(11, 3, 16)
    @test s9 == zeros(11, 3, 24)
    @test s10 == zeros(21, 4, 16)
    @test s3ma == repeat([_r1  _r2; 0  0; 0  0;;;]; outer=(1, 1, 10))
    @test s3mb == repeat([1  1; 1  1; 1  1;;;]; outer=(1, 1, 10))
    @test s3mc == repeat([1  1; 1  1; 1  1;;;]; outer=(1, 1, 3))
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
    @test_throws MethodError samplerenewaldidinfections(
        [0, 1], df3, 4.0:6;
        interventions=zeros(2, 2), 
        Ns=(100 .* ones(2)), 
        seedmatrix=[0  0; 1  1], 
        ngroups=2, 
        ntimes=2,
    )
    @test_throws MethodError samplerenewaldidinfections(
        [0, 1], df3, 4:0.5:6;
        interventions=zeros(2, 2), 
        Ns=(100 .* ones(2)), 
        seedmatrix=[0  0; 1  1], 
        ngroups=2, 
        ntimes=2,
    )
    @test_throws BoundsError samplerenewaldidinfections(
        [0, 1], df3, 4:16;  # df3 is 10 rows long
        interventions=zeros(2, 2), 
        Ns=(100 .* ones(2)), 
        seedmatrix=[0  0; 1  1], 
        ngroups=2, 
        ntimes=2,
    )
        @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1;
        interventions=zeros(10, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 4),
        ngroups=3,
        ntimes=10,
    )
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1;
        interventions=zeros(10, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(8, 3), 
        ngroups=3, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1;
        interventions=zeros(10, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(6, 3), 
        ngroups=3, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1;
        interventions=zeros(10, 2), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 3), 
        ngroups=2, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1;
        interventions=zeros(10, 4), 
        Ns=100 .* ones(4), 
        seedmatrix=zeros(7, 3), 
        ngroups=4, 
        ntimes=10,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1;
        interventions=zeros(9, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 3), 
        ngroups=3, 
        ntimes=9,
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(
        zeros(2), df1;
        interventions=zeros(11, 3), 
        Ns=100 .* ones(3), 
        seedmatrix=zeros(7, 3), 
        ngroups=3, 
        ntimes=11,
    ) 
end

@testset "quantiles of a single sample return a warning and the input" begin
    wm = "only a single sample provided to `quantilerenewaldidinfections`; input returned \
        unchanged for any value of `q`"
    @test quantilerenewaldidinfections(s1, 0.5; mutewarnings=true) == zeros(11, 3)
    @test_warn wm quantilerenewaldidinfections(s1, 0.5)
    @test_nowarn quantilerenewaldidinfections(s1, 0.5; mutewarnings=true)
    @test_warn wm quantilerenewaldidinfections(s1, 0.5; mutewarnings=false)
    @test quantilerenewaldidinfections(s2, 0.5; mutewarnings=true) == zeros(13, 4)
    @test quantilerenewaldidinfections(s5, 0.5; mutewarnings=true) == [1  1; 2  0.5; 4  0.25]
    @test quantilerenewaldidinfections(zeros(11, 3, 1), 0.5; mutewarnings=true) == zeros(11, 3)
    @test_warn wm quantilerenewaldidinfections(zeros(11, 3, 1), 0.5)
    @test_nowarn quantilerenewaldidinfections(zeros(11, 3, 1), 0.5; mutewarnings=true)
    @test_warn wm quantilerenewaldidinfections(zeros(11, 3, 1), 0.5; mutewarnings=false)
    @test quantilerenewaldidinfections(ones(20, 2, 1), 0.5; mutewarnings=true) == ones(20, 2)
    @test_warn wm quantilerenewaldidinfections(ones(20, 2, 1), 0.5)
    @test_nowarn quantilerenewaldidinfections(ones(20, 2, 1), 0.5; mutewarnings=true)
    @test_warn wm quantilerenewaldidinfections(ones(20, 2, 1), 0.5; mutewarnings=false)
    # no method currently intended for 4-dimensional arrays 
    @test_throws MethodError quantilerenewaldidinfections(ones(20, 2, 1, 1), 0.5)
    @test quantilerenewaldidinfections(s1, [0.25, 0.5, 0.75]; mutewarnings=true) == zeros(11, 3)
    @test_warn wm quantilerenewaldidinfections(s1, [0.25, 0.5, 0.75])
    @test_nowarn quantilerenewaldidinfections(s1, [0.25, 0.5, 0.75]; mutewarnings=true)
    @test_warn wm quantilerenewaldidinfections(s1, [0.25, 0.5, 0.75]; mutewarnings=false)
end

@testset "quantiles from arrays of samples" begin
    @test quantilerenewaldidinfections(s8, 0.5) == zeros(11, 3)
    @test_nowarn quantilerenewaldidinfections(s8, 0.5)
    @test_nowarn quantilerenewaldidinfections(s8, 0.5; mutewarnings=true)
    @test_nowarn quantilerenewaldidinfections(s8, 0.5; mutewarnings=false)
    A1 = rand(20, 3, 100)
    q025 = quantilerenewaldidinfections(A1, 0.025)
    q05 = quantilerenewaldidinfections(A1, 0.05)
    q25 = quantilerenewaldidinfections(A1, 0.25)
    q5 = quantilerenewaldidinfections(A1, 0.5)
    q75 = quantilerenewaldidinfections(A1, 0.75)
    q95 = quantilerenewaldidinfections(A1, 0.95)
    q975 = quantilerenewaldidinfections(A1, 0.975)
    @testset for t in 1:20, j in 1:3
        @test q025[t, j] < q05[t, j] < q25[t, j] < q5[t, j] < q75[t, j] < q95[t, j] < q975[t, j]
    end
    @test quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975]) == zeros(11, 3, 3)
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975])
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975]; mutewarnings=true)
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975]; mutewarnings=false)
    qv = quantilerenewaldidinfections(A1, [0.25, 0.5, 0.975])
    @test size(qv) == (20, 3, 3)
    @testset for (i, M) in enumerate([q25, q5, q975])
        @test qv[:, :, i] == M
    end
    @test_nowarn quantilerenewaldidinfections(A1, [0.25, 0.5, 0.975])
end
