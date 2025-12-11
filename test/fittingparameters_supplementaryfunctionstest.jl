# test functions that aid `RenewalDiD._renewaldid`

using Distributions: Normal
using Random: Xoshiro
using RenewalDiD
using StatsBase: mean, var
using Test

A1 = InterventionArray(4, [2, 3, nothing], Vector{Nothing}(nothing, 3))
A2 = InterventionArray(5, [2, 3], [4, nothing])
A3 = InterventionArray(4, [2, 3, nothing], [3, nothing, nothing])

obs1 = (_obs = zeros(Int, 20, 2); _obs[1, :] .+= 1; _obs[1:5, :] .+= 1; _obs)
obs2 = (_obs = zeros(Int, 20, 2); _obs[1, 1] += 1; _obs[1:5, :] .+= 1; _obs)

R0prediction1 = [
    0  1  2
    0  1  2
    0  1  2
    0  1  2
]
R0prediction2 = [
    -1    0    1
     0.5  1.5  2.5
     2.5  3.5  4.5
     0    1    2
]
R0prediction3 = [
    -1    0    1
    -0.5  1.5  2.5
     1.5  2.5  4.5
    -1    0    2
]
R0prediction4 = [0  0; 2  0; 2  2; 4  2; 4  2]
seedexpectation1 = [
    0.6       0.6
    0.689219  0.689219
    0.791705  0.791705
    0.90943   0.90943
    1.04466   1.04466
]
seedexpectation3 = [
    0.848528  0.848528
    0.90943   0.90943
    0.974703  0.974703
    1.04466   1.04466
    1.11964   1.11964
]
seedexpectation5 = [
    0.666667  0.666667
    0.765799  0.765799
    0.879672  0.879672
    1.01048   1.01048
    1.16073   1.16073
]
seedexpectation6 = [
    0.6       0.5
    0.689219  0.574349
    0.791705  0.659754
    0.90943   0.757858
    1.04466   0.870551
]
mseedexpectation01 = (m = zeros(7, 2); m[7, :] .+= 1; m)
mseedexpectation02 = (m = zeros(5, 2); m[5, :] .+= 2; m)
mseedexpectation03 = (m = zeros(5, 3); m[5, :] .+= 0.5; m)
mseedexpectation1 = [
    0.6       0.6
    0.689219  0.689219
    0.791705  0.791705
    0.90943   0.90943
    1.50965   1.50965
]
mseedexpectation2 = mseedexpectation1
mseedexpectation3 = seedexpectation3
mseedexpectation4 = [0  0; 0  0; 0  0; 0  0; 1  1]
mseedexpectation5 = (
    r4 = log(4 / 3) - 2 * log(2) / 5;
    r5 = 0.5 - r4;
    [0  0; 0  0; 0  0; r4  r4; r5  r5]
) 
mseedexpectation6 = [0  0; 0  0; 0  0; 0  0; 0.75  0.75]
predictedinfections3 = [1  0; 2  1; 0  0; 0  0; 0  0]
predictedinfections4 = [1  0; 2  1; 3  1.5; 4.5  2.25; 6.75  3.375]
predictedinfections5 = [3  0; 1  2; 2  3; 0  6; 0  12]
predictedinfections8 = ComplexF64[
    1+0.9im  0+1im; 
    2+0.7im  1+0.9im; 
    0+0.7im  0+0.9im; 
    0+0.7im  0+0.9im; 
    0+0.7im  0+0.9im
]
predictedinfections9 = let 
    inf3_1 = 3 * 0.7 
    sus3_1 = 0.7 - inf3_1 / 10 
    inf4_1 = 1.5 * inf3_1 * sus3_1 
    sus4_1 = sus3_1 - inf4_1 / 10 
    inf5_1 = 1.5 * inf4_1 * sus4_1 
    sus5_1 = sus4_1 - inf5_1 / 10 
    inf3_2 = 1.5 * 0.9 
    sus3_2 = 0.9 - inf3_2 / 10 
    inf4_2 = 1.5 * inf3_2 * sus3_2 
    sus4_2 = sus3_2 - inf4_2 / 10 
    inf5_2 = 1.5 * inf4_2 * sus4_2 
    sus5_2 = sus4_2 - inf5_2 / 10 
    ComplexF64[
        1+0.9im  0+1im; 
        2+0.7im  1+0.9im; 
        inf3_1+sus3_1*im  inf3_2+sus3_2*im; 
        inf4_1+sus4_1*im  inf4_2+sus4_2*im; 
        inf5_1+sus5_1*im  inf5_2+sus5_2*im
    ]
end

predictedcalcseedinfections4 = [3  0; 1  2]
predictedcalcseedinfections8 = ComplexF64[1+0.9im  0+1im; 2+0.7im  1+0.9im]
predictedcalcseedinfections9 = ComplexF64[1+0.98im  0+1im; 2+0.94im  1+0.99im]
predictedcalcseedinfections10 = ComplexF64[3+0.94im  0+1im; 1+0.92im  2+0.98im]

zerosstruct = RenewalDiDData( ; 
    observedcases=zeros(4, 2), interventions=zeros(3, 2), Ns=zeros(Int, 2)
)
namedzerosstruct = RenewalDiDData( ; 
    observedcases=zeros(4, 2), interventions=zeros(3, 2), Ns=zeros(Int, 2), id="IHaveAName"
)
namedzerostructunlimitied = RenewalDiDData( ; 
    observedcases=zeros(4, 2), 
    interventions=zeros(3, 2), 
    id="InfinitePopulation"
)

nzoutput = "RenewalDiDData{Float64, Matrix{Float64}, Vector{Int64}}\n observedcases:  [0.0 0.0; 0.0 0.0; \
    0.0 0.0; 0.0 0.0]\n interventions:  [0.0 0.0; 0.0 0.0; 0.0 0.0]\n Ns:             \
    [0, 0]\n exptdseedcases: [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.5 0.5]"
nnzoutput = "RenewalDiDData{Float64, Matrix{Float64}, Vector{Int64}}, (IHaveAName)\n observedcases:  [0.0 \
    0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]\n interventions:  [0.0 0.0; 0.0 0.0; 0.0 0.0]\n \
    Ns:             [0, 0]\n exptdseedcases: [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; \
    0.0 0.0; 0.5 0.5]"
unnzoutput = "RenewalDiDData{Float64, Matrix{Float64}, Nothing}, \
    (InfinitePopulation)\n observedcases:  [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 \
    0.0]\n interventions:  [0.0 0.0; 0.0 0.0; 0.0 0.0]\n Ns:             \
    unlimited\n exptdseedcases: [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.5 0.5]"

@testset "data struct" begin
    @test (@inferred RenewalDiDData( ; 
        observedcases=zeros(4, 2), interventions=zeros(3, 2), Ns=zeros(Int, 2)
    )) == zerosstruct
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(3, 2), interventions=zeros(3, 2), Ns=zeros(Int, 2)
    )
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(4, 3), interventions=zeros(3, 2), Ns=zeros(Int, 2)
    )
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(4, 2), interventions=zeros(3, 3), Ns=zeros(Int, 2)
    )
    @test_throws DimensionMismatch RenewalDiDData( ; 
        observedcases=zeros(4, 2), interventions=zeros(3, 2), Ns=zeros(Int, 3)
    )
end

@testset "calculate ntimes" begin
    @test RenewalDiD._ntimes(zeros(2)) == 2
    @test RenewalDiD._ntimes(zeros(2, 3)) == 2
    @test RenewalDiD._ntimes(zeros(3, 2)) == 3
end

@testset "calculate ngroups" begin
    @test RenewalDiD._ngroups(zeros(2)) == 1
    @test RenewalDiD._ngroups(zeros(2, 3)) == 3
    @test RenewalDiD._ngroups(zeros(3, 2)) == 2
end

@testset "calculate ninterventions" begin
    @test RenewalDiD._ninterventions(zeros(2, 3)) == 1
    @test RenewalDiD._ninterventions(zeros(3, 2)) == 1
    @test RenewalDiD._ninterventions(zeros(3, 2, 6)) == 6
    @test RenewalDiD._ninterventions(A1) == 2
end

@testset "generate vector of gamma values" begin
    @test RenewalDiD._gammavec(zeros(3), 1) == zeros(4)
    @test RenewalDiD._gammavec(zeros(4), 1) == zeros(5)
    @test RenewalDiD._gammavec([1, -1, 0], 1) == [1, -1, 0, 0]
    @test RenewalDiD._gammavec([1, -1, 0], 0.5) == [0.5, -0.5, 0, 0]
    @test RenewalDiD._gammavec([1, -2, 0.5], 0.5) == [0.5, -1, 0.25, 0.25]
    @test RenewalDiD._gammavec([3, -2, 0.5], 0.5) == [1.5, -1, 0.25, -0.75]
    g = RenewalDiD._gammavec(rand(3), rand())
    @test length(g) == 4
    @test sum(g) ≈ 0 atol=1e-10
end

@testset "expected variance of gamma values" begin
    rng = Xoshiro(1)
    gammasigma = 0.5 
    vecofgammavec = [
        RenewalDiD._gammavec(rand(rng, Normal(0, 1), 3), gammasigma) 
        for _ in 1:10_000
    ]
    vecofgammavars = [var(vecofgammavec[i]) for i in 1:10_000]
    @test mean(vecofgammavars) ≈ gammasigma atol=1e-2
end

@testset "generate vector of theta values" begin
    @test RenewalDiD._assemblethetavec(zeros(3), 1, 4, RenewalDiD.automatic) == zeros(4)
    @test RenewalDiD._assemblethetavec(zeros(3), 1, 4, 1) == zeros(4)
    @test RenewalDiD._assemblethetavec(zeros(3), 1, 4, 1, 0) == zeros(4)
    @test_throws DimensionMismatch RenewalDiD._assemblethetavec(zeros(3), 1, 3, 1)
    @test_throws ArgumentError RenewalDiD._assemblethetavec(zeros(3), 1, 4, 0)
    @test_throws ArgumentError RenewalDiD._assemblethetavec(zeros(3), 1, 4, -1)
    @test RenewalDiD._assemblethetavec([1, 2, 3], 1, 6, 2) == [0, 0.5, 1, 1.5, 2, 3] 
    @test RenewalDiD._assemblethetavec([1, 2, 3], 1, 7, 2) == 0:0.5:3
    @test_throws DimensionMismatch RenewalDiD._thetavec(zeros(3), 1, 3)
    @test RenewalDiD._thetavec(zeros(3), 1, 4) == zeros(4)
    @test RenewalDiD._thetavec(zeros(4), 1, 5) == zeros(5)
    @test RenewalDiD._thetavec(ones(3), 1, 4) == [0, 1, 2, 3]
    @test RenewalDiD._thetavec(ones(3), 0.5, 4) == [0, 0.5, 1, 1.5]
    @test RenewalDiD._thetavec([-2.5, 0.5, 0.5], 0.5, 4) == [0, -1.25, -1, -0.75]
    @test RenewalDiD._thetavec([-2.5, 0.5, 0.5], 0.5, 4; thetainterval=RenewalDiD.automatic) == 
        [0, -1.25, -1, -0.75]
    @test RenewalDiD._thetavec([-2.5, 0.5, 0.5], 0.5, 4; thetainterval=1) == 
        [0, -1.25, -1, -0.75]
    @test RenewalDiD._thetavec([1, 2, 3], 1, 6; thetainterval=2) == 
        [0, 0.5, 1.5, 3, 5, 8] 
    @test RenewalDiD._thetavec([1, 2, 3], 1, 7; thetainterval=2) == 
        [0, 0.5, 1.5, 3, 5, 7.5, 10.5]
    @test RenewalDiD._thetavec([1, 2, 3], 0.5, 7; thetainterval=2) == 
        [0, 0.25, 0.75, 1.5, 2.5, 3.75, 5.25]
    # test length as shorthand for testing it successfully generates vector
    @test length(RenewalDiD._thetavec(rand(15), rand(), 100; thetainterval=7)) == 100
    @test_throws ArgumentError RenewalDiD._thetavec(rand(14), rand(), 100; thetainterval=7)
    @test_throws ArgumentError RenewalDiD._thetavec(rand(16), rand(), 100; thetainterval=7)
end

@testset "generate matrix of log R_0 values with multiple interventions" begin
    @test RenewalDiD._predictedlogR_0(0, zeros(3), zeros(4), zeros(2), A1) == zeros(4, 3)
    @test RenewalDiD._predictedlogR_0(0, zeros(2), zeros(5), zeros(2), A2) == zeros(5, 2)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(2), zeros(4), 0, A1)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(3), zeros(3), 0, A1)
    @test RenewalDiD._predictedlogR_0(1, zeros(3), zeros(4), zeros(2), A1) == ones(4, 3)
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], zeros(4), zeros(2), A1) == R0prediction1
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], zeros(2), A1) == R0prediction2
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], [-1, -1], A1) == R0prediction3
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], zeros(4), zeros(2), A3) == R0prediction1
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], zeros(2), A3) == R0prediction2
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], [-1, 0], A3) == R0prediction3
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], [-1, -1], A3) != R0prediction3
    @test RenewalDiD._predictedlogR_0(0, zeros(2), zeros(5), [2, 2], A2) == R0prediction4
end

@testset "output from `Base.show`" begin
    @test repr("text/plain", zerosstruct) == nzoutput 
    @test repr("text/plain", namedzerosstruct) == nnzoutput 
    @test repr("text/plain", namedzerostructunlimitied) == unnzoutput 
end
