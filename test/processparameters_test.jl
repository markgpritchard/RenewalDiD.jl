# test functions for processing output of fitted parameters

# `minsigma2 = 0` for all tests - may wish to revisit these tests, or even how `minsigma2` 
# is used

using RenewalDiD
using RenewalDiD: testdataframe, testsimulation
using RenewalDiD: tupleforsamplerenewaldidinfections as tfsr
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
    taudefault=[0, rand()], 
    alpha=[0, rand()], 
    sigma_gamma=[0, rand()], 
    sigma_theta=[0, rand()], 
    minsigma2=zeros(2),
)
df2 = testdataframe( ; 
    nchains=4, 
    niterations=200, 
    ngroups=4, 
    ntimes=12, 
    nseeds=7,
    minsigma2=zeros(800),
)
df3 = testdataframe( ;
    nchains=1,
    niterations=10,
    ngroups=2,
    ntimes=2,
    nseeds=2,
    alpha=zeros(10),
    gammadefault=zeros(10), 
    psi=ones(10),
    thetadefault=zeros(10), 
    mxdefault=zeros(10),
    minsigma2=zeros(10),
)
df4 = testdataframe( ;
    nchains=1,
    niterations=1,
    ngroups=2,
    ntimes=2,
    nseeds=2,
    alpha=[log(2)],
    gammadefault=zeros(1), 
    psi=ones(1),
    thetadefault=zeros(1), 
    mxdefault=zeros(1),
    minsigma2=zeros(1),
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
        psi=ones(1),
        thetadefault=zeros(1), 
        mxdefault=zeros(1),
        minsigma2=zeros(1),
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
        psi=ones(1),
        thetadefault=zeros(1), 
        mxdefault=zeros(1),
        minsigma2=zeros(1),
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
        minsigma2=zeros(1),
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
    minsigma2=zeros(16),
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
    minsigma2=zeros(24),
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
    minsigma2=zeros(16),
)
data1 = RenewalDiDData( ; 
    observedcases=zeros(11, 3), 
    interventions=zeros(10, 3), 
    exptdseedcases=zeros(7, 3),
)
data2 = RenewalDiDData( ; 
    observedcases=zeros(13, 4), 
    interventions=zeros(12, 4), 
    exptdseedcases=zeros(7, 4), 
)
data3 = RenewalDiDData( ; 
    observedcases=zeros(3, 2), 
    interventions=zeros(2, 2), 
    exptdseedcases=[0  0; 1  1],
)
data8 = RenewalDiDData( ; 
    observedcases=zeros(11, 3), 
    interventions=zeros(10, 3), 
    exptdseedcases=zeros(2, 3),
)
data10 = RenewalDiDData( ; 
    observedcases=zeros(21, 4), 
    interventions=zeros(20, 4), 
    exptdseedcases=zeros(2, 4),
)
data3ma = RenewalDiDData( ; 
    observedcases=zeros(3, 2), 
    interventions=zeros(2, 2), 
    exptdseedcases=[0  0; _r1  _r2],
)
data3mb = RenewalDiDData( ; 
    observedcases=zeros(3, 2), 
    interventions=zeros(2, 2), 
    exptdseedcases=[0  0; 1  1],
)
s1 = samplerenewaldidinfections(tfsr(data1; vec=zeros(2)), df1, 2).output
s2 = samplerenewaldidinfections(tfsr(data2; vec=zeros(2)), df2, 2).output
s3 = samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df3, 1).output
s4 = samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df4, 1).output
s5 = samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df5, 1).output
s6 = samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df6, 1).output
s7 = samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df7, 1).output
s8 = samplerenewaldidinfections(tfsr(data8; vec=zeros(2)), df8).output
s9 = samplerenewaldidinfections(tfsr(data8; vec=zeros(2)), df9).output
s10 = samplerenewaldidinfections(tfsr(data10; vec=zeros(2)), df10).output
s3ma = samplerenewaldidinfections(tfsr(data3ma; vec=zeros(2)), df3).output
s3mb = samplerenewaldidinfections(tfsr(data3mb; vec=[0, 1]), df3).output
s3mc = samplerenewaldidinfections(tfsr(data3mb; vec=[0, 1]), df3, 4:6).output
s1r1 = samplerenewaldidinfections(tfsr(data1; vec=zeros(2)), df1, 2; repeatsamples=2).output
s1r2 = samplerenewaldidinfections(tfsr(data1; vec=zeros(2)), df1, 2; repeatsamples=5).output
s2r1 = samplerenewaldidinfections(tfsr(data2; vec=zeros(2)), df2, 2; repeatsamples=2).output
s2r2 = samplerenewaldidinfections(tfsr(data2; vec=zeros(2)), df2, 2; repeatsamples=10).output
s3r = samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df3, 1; repeatsamples=4).output
s8r = samplerenewaldidinfections(tfsr(data8; vec=zeros(2)), df8; repeatsamples=3).output
s3mar = samplerenewaldidinfections(tfsr(data3ma; vec=zeros(2)), df3; repeatsamples=3).output
s3mbr = samplerenewaldidinfections(tfsr(data3mb; vec=[0, 1]), df3).output
s3mcr = samplerenewaldidinfections(tfsr(data3mb; vec=[0, 1]), df3, 4:6).output

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

@testset "sample multiple rows" begin
    @test s8 == zeros(11, 3, 16)
    @test s9 == zeros(11, 3, 24)
    @test s10 == zeros(21, 4, 16)
    @test s3ma == repeat([_r1  _r2; 0  0; 0  0;;;]; outer=(1, 1, 10))
    @test s3mb == repeat([1  1; 1  1; 1  1;;;]; outer=(1, 1, 10))
    @test s3mc == repeat([1  1; 1  1; 1  1;;;]; outer=(1, 1, 3))
end

@testset "sampling errors" begin
    # tests updated as some dimension mismatch errors should now be caught while calling 
    # `RenewalDiDData`
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(11, 3), 
        interventions=zeros(10, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 4),
    )
    # tests removed that used to throw an error for an older form of this function 
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(11, 2), 
        interventions=zeros(10, 2), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 3),
    )
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(11, 4), 
        interventions=zeros(10, 4), 
        Ns=(100 .* ones(4)), 
        exptdseedcases=zeros(7, 3),
    )
    data13 = RenewalDiDData( ; 
        observedcases=zeros(10, 3), 
        interventions=zeros(9, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 3),
    )
    @test_throws DimensionMismatch samplerenewaldidinfections(tfsr(data13; vec=zeros(2)), df1, 2) 
    data14 = RenewalDiDData( ; 
        observedcases=zeros(12, 3), 
        interventions=zeros(11, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 3),
    )
    @test_throws DimensionMismatch samplerenewaldidinfections(tfsr(data14; vec=zeros(2)), df1, 2) 
    @test_throws BoundsError samplerenewaldidinfections(tfsr(data1; vec=zeros(2)), df1, 3) 
    @test_throws MethodError samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df3, 4.0:6)
    @test_throws MethodError samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df3, 4:0.5:6)
    @test_throws BoundsError samplerenewaldidinfections(tfsr(data3; vec=[0, 1]), df3, 4:16)  # df3 is 10 rows long
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(11, 3), 
        interventions=zeros(10, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 4),
    )
    data16 = RenewalDiDData( ; 
        observedcases=zeros(11, 3), 
        interventions=zeros(10, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(8, 3),
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(tfsr(data16; vec=zeros(2)), df1) 
    data17 = RenewalDiDData( ; 
        observedcases=zeros(11, 3), 
        interventions=zeros(10, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(6, 3),
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(tfsr(data17; vec=zeros(2)), df1) 
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(11, 2), 
        interventions=zeros(10, 2), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 3),
    ) 
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(11, 4), 
        interventions=zeros(10, 4), 
        Ns=(100 .* ones(4)), 
        exptdseedcases=zeros(7, 3),
    ) 
    data18 = RenewalDiDData( ; 
        observedcases=zeros(10, 3), 
        interventions=zeros(9, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 3),
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(tfsr(data18; vec=zeros(2)), df1) 
    data19 = RenewalDiDData( ; 
        observedcases=zeros(12, 3), 
        interventions=zeros(11, 3), 
        Ns=(100 .* ones(3)), 
        exptdseedcases=zeros(7, 3),
    ) 
    @test_throws DimensionMismatch samplerenewaldidinfections(tfsr(data19; vec=zeros(2)), df1) 
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
    wm1 = "(0.25, 0.975): other functions expect that credible intervals are symmetrical"
    wm2 = "0.05 < 0.25: other functions expect quantiles are calculated in ascending order"
    wm3 = "[0.25, 0.5, 0.75, 0.975]: other functions expect an odd number of quantiles"
    wm4 = "[0.25, 0.75, 0.975]: other functions expect the middle quantile is the median"
    wm5 = "[0.25, 0.75, 0.5]: other functions expect the middle quantile is the median"
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
    @test quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975]; mutewarnings=true) == zeros(11, 3, 3)
    @test_warn wm1 quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975])
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975]; mutewarnings=true)
    @test_warn wm1 quantilerenewaldidinfections(s8, [0.25, 0.5, 0.975]; mutewarnings=false)
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.75])
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.75]; mutewarnings=true)
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.75]; mutewarnings=false)
    @test_nowarn quantilerenewaldidinfections(s8, [0.05, 0.25, 0.5, 0.75, 0.95])
    @test_nowarn quantilerenewaldidinfections(s8, [0.05, 0.25, 0.5, 0.75, 0.95]; mutewarnings=true)
    @test_nowarn quantilerenewaldidinfections(s8, [0.05, 0.25, 0.5, 0.75, 0.95]; mutewarnings=false)
    @test_warn wm2 quantilerenewaldidinfections(s8, [0.25, 0.05, 0.5, 0.95, 0.75])
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.05, 0.5, 0.95, 0.75]; mutewarnings=true)
    @test_warn wm2 quantilerenewaldidinfections(s8, [0.25, 0.05, 0.5, 0.95, 0.75]; mutewarnings=false)
    @test_warn wm3 quantilerenewaldidinfections(s8, [0.25, 0.5, 0.75, 0.975])
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.5, 0.75, 0.975]; mutewarnings=true)
    @test_warn wm3 quantilerenewaldidinfections(s8, [0.25, 0.5, 0.75, 0.975]; mutewarnings=false)
    @test_warn wm4 quantilerenewaldidinfections(s8, [0.25, 0.75, 0.975])
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.75, 0.975]; mutewarnings=true)
    @test_warn wm4 quantilerenewaldidinfections(s8, [0.25, 0.75, 0.975]; mutewarnings=false)
    @test_warn wm5 quantilerenewaldidinfections(s8, [0.25, 0.75, 0.5])
    @test_nowarn quantilerenewaldidinfections(s8, [0.25, 0.75, 0.5]; mutewarnings=true)
    @test_warn wm5 quantilerenewaldidinfections(s8, [0.25, 0.75, 0.5]; mutewarnings=false)
    qv = quantilerenewaldidinfections(A1, [0.25, 0.5, 0.975]; mutewarnings=true)
    @test size(qv) == (20, 3, 3)
    @testset for (i, M) in enumerate([q25, q5, q975])
        @test qv[:, :, i] == M
    end
    @test_nowarn quantilerenewaldidinfections(A1, [0.25, 0.5, 0.75])
end

@testset "repeat samples" begin
    @test samplerenewaldidinfections(tfsr(data1; vec=zeros(2)), df1, 2; repeatsamples=nothing).output == s1
    @testset for i in eachindex(s1r1)
        @test s1r1[i] ≈ 0 atol=1e-3 
    end
    @test s1r1 == zeros(11, 3, 2)
    @test s1r2 == zeros(11, 3, 5)
    @test s2r1 == zeros(13, 4, 2)
    @test s2r2 == zeros(13, 4, 10)
    @test s3r == ones(3, 2, 4)
    @test s8r == zeros(11, 3, 48)
    @test s3mar == repeat([_r1  _r2; 0  0; 0  0;;;]; outer=(1, 1, 30))
end

@testset "intervention starttimes" begin
    _im1 = InterventionMatrix(100, [40, 80, nothing, -1, 101])
    _im2 = InterventionMatrix{Bool}(100, [40, 80, nothing, -1, 101])
    _am1 = let 
        m = zeros(100, 5)
        for t in 40:100 
            m[t, 1] = 1
            t < 80 && continue 
            m[t, 2] = 1
        end 
        m
    end
    _am2 = let 
        m = zeros(Bool, 100, 5)
        for t in 40:100 
            m[t, 1] = true
            t < 80 && continue 
            m[t, 2] = true
        end 
        m
    end
    @testset for M in [_im1, _im2, _am1, _am2]
        @test RenewalDiD._interventionstarttimes(M, 1) == 40 
        @test RenewalDiD._interventionstarttimes(M, 2) == 80
        @testset for i in 3:5 
            @test isnothing(RenewalDiD._interventionstarttimes(M, i))
        end
        @test RenewalDiD._interventionstarttimes(M) == [40, 80, nothing, nothing, nothing]
    end
end
