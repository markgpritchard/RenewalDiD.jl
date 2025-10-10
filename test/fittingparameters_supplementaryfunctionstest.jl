# test functions that aid `RenewalDiD._renewaldid`

using RenewalDiD
using Test
using Random
using StatsBase: mean, var
using Turing: Normal

M1 = InterventionMatrix(4, [2, 3, nothing]) 
M2 = InterventionMatrix(5, Vector{Nothing}(nothing, 2); mutewarnings=true) 
A1 = InterventionArray(4, [2, 3, nothing], Vector{Nothing}(nothing, 3))
A2 = InterventionArray(5, [2, 3], [4, nothing])
A3 = InterventionArray(4, [2, 3, nothing], [3, nothing, nothing])

obs1 = (_obs = zeros(Int, 20, 2); _obs[1, :] .+= 1; _obs[1:5, :] .+= 1; _obs)
obs2 = (_obs = zeros(Int, 20, 2); _obs[1, 1] += 1; _obs[1:5, :] .+= 1; _obs)

seedinfections1 = [1  0; 2  1]
M_x1 = [2  0; -0.5  1; 1  0.5; -2  1; 0  1]

calcseedinfections1 = let  # note, Ns not yet used  
    infectionsmatrix = zeros(2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(5, 2), zeros(2, 2), nothing, 2)
    infectionsmatrix
end
calcseedinfections2 = let  # note, Ns not yet used  
    infectionsmatrix = zeros(3, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(6, 2), zeros(3, 2), nothing, 3)
    infectionsmatrix
end
calcseedinfections3 = let  # note, Ns not yet used  
    infectionsmatrix = zeros(2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(5, 2), seedinfections1, nothing, 2)
    infectionsmatrix
end
calcseedinfections4 = let  # note, Ns not yet used  
    infectionsmatrix = zeros(2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, M_x1, seedinfections1, nothing, 2)
    infectionsmatrix
end
calcseedinfections5 = let 
    infectionsmatrix = zeros(ComplexF64, 2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(5, 2), zeros(2, 2), ones(2), 2)
    infectionsmatrix
end
calcseedinfections6 = let 
    infectionsmatrix = zeros(ComplexF64, 2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(5, 2), zeros(2, 2), [1000, 2000], 2)
    infectionsmatrix
end
calcseedinfections7 = let 
    infectionsmatrix = zeros(ComplexF64, 3, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(6, 2), zeros(3, 2), ones(2), 3)
    infectionsmatrix
end
calcseedinfections8 = let 
    infectionsmatrix = zeros(ComplexF64, 2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(5, 2), seedinfections1, [10, 10], 2)
    infectionsmatrix
end
calcseedinfections9 = let 
    infectionsmatrix = zeros(ComplexF64, 2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, zeros(5, 2), seedinfections1, [50, 100], 2)
    infectionsmatrix
end
calcseedinfections10 = let  
    infectionsmatrix = zeros(ComplexF64, 2, 2)
    RenewalDiD._infections_seed!(infectionsmatrix, M_x1, seedinfections1, [50, 100], 2)
    infectionsmatrix
end
calcinfections1 = RenewalDiD._infections(
    g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), nothing, 2
)
calcinfections2 = RenewalDiD._infections(
    g_covid, ComplexF64, zeros(6, 2), log.(zeros(4, 2)), zeros(2, 2), nothing, 2
)
calcinfections3 = RenewalDiD._infections(
    generationtime, 
    ComplexF64, 
    zeros(5, 2), 
    log.(zeros(3, 2)), 
    seedinfections1,  
    nothing,  # tests below assume no change in susceptibility 
    2;
    vec=[0, 1],
)
calcinfections4 = RenewalDiD._infections(
    generationtime, 
    ComplexF64, 
    zeros(5, 2), 
    log.(1.5 .* ones(3, 2)), 
    seedinfections1, 
    nothing,  # tests below assume no change in susceptibility 
    2;
    vec=[0, 1],
)
calcinfections5 = RenewalDiD._infections(
    generationtime, 
    ComplexF64, 
    M_x1, 
    log.(ones(3, 2)), 
    seedinfections1, 
    nothing,  # tests below assume no change in susceptibility 
    2;
    vec=[0, 1]
)
calcinfections6 = RenewalDiD._infections(
    g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), [10, 10], 2
)
calcinfections7 = RenewalDiD._infections(
    g_covid, ComplexF64, zeros(6, 2), log.(zeros(4, 2)), zeros(2, 2), [10, 10], 2
)
calcinfections8 = RenewalDiD._infections(
    generationtime, 
    ComplexF64, 
    zeros(5, 2), 
    log.(zeros(3, 2)), 
    seedinfections1, 
    [10, 10], 
    2;
    vec=[0, 1],
)
calcinfections9 = RenewalDiD._infections(
    generationtime, 
    ComplexF64, 
    zeros(5, 2), 
    log.(1.5 * ones(3, 2)), 
    seedinfections1, 
    [10, 10], 
    2;
    vec=[0, 1],
)

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
    0  0; 0  0; 0  0; 0  0; (log(6 / 5) - log(2) / 5)  (log(6 / 5) - log(2) / 5)
]
seedexpectation2 = [
    0  0; 0  0; 0  0; 0  0; (log(6 / 5) - log(2) / 5)  (log(6 / 5) - log(2) / 5)
]
seedexpectation3 = [
    0  0; 
    0  0; 
    0  0; 
    (log(6 / 5) - 2 * log(2) / 10)  (log(6 / 5) - 2 * log(2) / 10); 
    (log(6 / 5) - log(2) / 10)  (log(6 / 5) - log(2) / 10)
]
seedexpectation4 = [
    0  0; 0  0; 0  0; 0  0; (log(6 / 5) - log(2) / 5)  (log(6 / 5) - log(2) / 5)
]
seedexpectation5 = [
    0  0; 
    0  0; 
    0  0; 
    (log(4 / 3) - 2 * log(2) / 5)  (log(4 / 3) - 2 * log(2) / 5); 
    (log(4 / 3) - log(2) / 5)  (log(4 / 3) - log(2) / 5)
]
seedexpectation6 = [0  0; 0  0; 0  0; 0  0; (log(6 / 5) - log(2) / 5)  0]
mseedexpectation01 = (m = zeros(7, 2); m[7, :] .+= 1; m)
mseedexpectation02 = (m = zeros(5, 2); m[5, :] .+= 2; m)
mseedexpectation03 = (m = zeros(5, 3); m[5, :] .+= 0.5; m)
mseedexpectation1 = [0  0; 0  0; 0  0; 0  0; 0.5  0.5]
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
    @test RenewalDiD._ntimes(M1) == 4
end

@testset "calculate ngroups" begin
    @test RenewalDiD._ngroups(zeros(2)) == 1
    @test RenewalDiD._ngroups(zeros(2, 3)) == 3
    @test RenewalDiD._ngroups(zeros(3, 2)) == 2
    @test RenewalDiD._ngroups(M1) == 3
end

@testset "calculate ninterventions" begin
    @test RenewalDiD._ninterventions(zeros(2, 3)) == 1
    @test RenewalDiD._ninterventions(zeros(3, 2)) == 1
    @test RenewalDiD._ninterventions(M1) == 1
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
    @test length(enewalDiD._thetavec(rand(15), rand(), 100; thetainterval=7)) == 100
end

@testset "generate matrix of log R_0 values" begin
    @test RenewalDiD._predictedlogR_0(0, zeros(3), zeros(4), [0], M1) == zeros(4, 3)
    @test RenewalDiD._predictedlogR_0(0, zeros(2), zeros(5), [2], M2) == zeros(5, 2)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(2), zeros(4), [0], M1)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(3), zeros(3), [0], M1)
    @test RenewalDiD._predictedlogR_0(1, zeros(3), zeros(4), [0], M1) == ones(4, 3)
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], zeros(4), [0], M1) == R0prediction1
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], [0], M1) == R0prediction2
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], [-1], M1) == R0prediction3
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

@testset "seed infections returned by `_infections_seed" begin
    @test RenewalDiD._infections_seed!(zeros(2, 2), zeros(5, 2), zeros(2, 2), nothing, 2) == 
        calcseedinfections1
    @test RenewalDiD._infections_seed!(
        zeros(ComplexF64, 2, 2), M_x1, seedinfections1, [50, 100], 2
    ) == calcseedinfections10
end

@testset "expected seed cases" begin
    @test expectedseedcases(zeros(20, 2), 7; minvalue=0) == zeros(7, 2)
    @test expectedseedcases(zeros(20, 2), 5; minvalue=0) == zeros(5, 2)
    @test expectedseedcases(zeros(20, 3), 5; minvalue=0) == zeros(5, 3)
    @test expectedseedcases(obs1, 5; minvalue=0) == seedexpectation1
    @test expectedseedcases(obs1, 5; doubletime=5, minvalue=0) == seedexpectation2
    @test expectedseedcases(obs1, 5; doubletime=10, minvalue=0) == seedexpectation3
    @test expectedseedcases(obs1, 5; sampletime=5, minvalue=0) == seedexpectation4
    @test expectedseedcases(obs1, 5; sampletime=3, minvalue=0) == seedexpectation5
    @test expectedseedcases(obs2, 5; minvalue=0) == seedexpectation6
    @test expectedseedcases(zeros(2, 3), 5; minvalue=0) == zeros(5, 3)
    @test_throws BoundsError expectedseedcases(zeros(2, 3), 5; sampletime=5)
    mm01 = expectedseedcases(zeros(20, 2), 7; minvalue=1)
    @testset for i in eachindex(mm01)
        @test mm01[i] ≈ mseedexpectation01[i] atol=1e-10
    end
    mm02 = expectedseedcases(zeros(20, 2), 5; minvalue=2)
    @testset for i in eachindex(mm02)
        @test mm02[i] ≈ mseedexpectation02[i] atol=1e-10
    end
    mm03 = expectedseedcases(zeros(20, 3), 5; minvalue=0.5) 
    @testset for i in eachindex(mm03)
        @test mm03[i] ≈ mseedexpectation03[i] atol=1e-10
    end
    mm03a = expectedseedcases(zeros(20, 3), 5)  # `minvalue=0.5` is now default
    @testset for i in eachindex(mm03)
        @test mm03a[i] ≈ mseedexpectation03[i] atol=1e-10
    end
    mm1 = expectedseedcases(obs1, 5; minvalue=0.5) 
    @testset for i in eachindex(mm1)
        @test mm1[i] ≈ mseedexpectation1[i] atol=1e-10
    end
    mm1a = expectedseedcases(obs1, 5) 
    @testset for i in eachindex(mm1)
        @test mm1a[i] ≈ mseedexpectation1[i] atol=1e-10
    end
    mm2 = expectedseedcases(obs1, 5; doubletime=5, minvalue=0.5) 
    @testset for i in eachindex(mm2)
        @test mm2[i] ≈ mseedexpectation2[i] atol=1e-10
    end
    mm2a = expectedseedcases(obs1, 5; doubletime=5) 
    @testset for i in eachindex(mm2)
        @test mm2a[i] ≈ mseedexpectation2[i] atol=1e-10
    end
    mm3 = expectedseedcases(obs1, 5; doubletime=10, minvalue=0.1) 
    @testset for i in eachindex(mm3)
        @test mm3[i] ≈ mseedexpectation3[i] atol=1e-10
    end
    mm4 = expectedseedcases(obs1, 5; sampletime=5, minvalue=1) 
    @testset for i in eachindex(mm4)
        @test mm4[i] ≈ mseedexpectation4[i] atol=1e-10
    end
    mm5 = expectedseedcases(obs1, 5; sampletime=3, minvalue=0.5) 
    @testset for i in eachindex(mm5)
        @test mm5[i] ≈ mseedexpectation5[i] atol=1e-10
    end
    mm5a = expectedseedcases(obs1, 5; sampletime=3) 
    @testset for i in eachindex(mm5)
        @test mm5a[i] ≈ mseedexpectation5[i] atol=1e-10
    end
    mm6 = expectedseedcases(obs2, 5; minvalue=0.75)
    @testset for i in eachindex(mm6)
        @test mm6[i] ≈ mseedexpectation6[i] atol=1e-10
    end
end
        
@testset "expected number of infections" begin
    # using `log(0)` and `log(1)` to emphasize that these parameters are natural logarithms
    @test RenewalDiD._expectedinfections(g_covid, log(0), zeros(10)) == 0
    @test RenewalDiD._expectedinfections(g_covid, log(1), ones(10)) == 
        sum(RenewalDiD.COVIDSERIALINTERVAL[2:11])
    @test RenewalDiD._expectedinfections(g_covid, log(1), ones(6)) == 
    sum(RenewalDiD.COVIDSERIALINTERVAL[2:7])
    @test RenewalDiD._expectedinfections(g_covid, log(1), [0, 0, 0, 0, 0, 1]) == 0.0440204506
    @test RenewalDiD._expectedinfections(g_covid, log(1), [1, 0, 0, 0, 0, 0]) == 0.0917470443
    @test RenewalDiD._expectedinfections(g_seir, log(0), zeros(10); mu=0.5) == 0
    @test RenewalDiD._expectedinfections(g_seir, log(1), [1, 0, 0, 0, 0, 0]; mu=0.5) == 
        0.07468060255179591
    @test_throws UndefKeywordError RenewalDiD._expectedinfections(g_seir, log(0), zeros(10))
    # we do not need to use `generationtime` with `g_seir` but test it as though `g_seir` 
    # were a user-generated function
    @test RenewalDiD._expectedinfections(
        generationtime, log(0), zeros(10); 
        func=g_seir, mu=0.5,
    ) == 0
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), [1, 0, 0, 0, 0, 0]; 
        func=g_seir, mu=0.5,
    ) == 0.07468060255179591
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), [1, 0, 0, 0, 0, 0]; 
        func=g_seir, mu=0.5, t_max=5,
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
        func=g_seir, mu=0.5, vec=[0, 0, 0, 0, 0, 1],
    ) 
end   
        
@testset "approximate number of cases" begin
    @testset "with no variance is just the point estimate" begin
        @test RenewalDiD._approxcases(0, 0) == 0
        @test RenewalDiD._approxcases(1, 0) == 1
    end
    @testset "account for variance" begin
        @test RenewalDiD._approxcases(1, 0.5) == 1.5
        @test RenewalDiD._approxcases(2, 0.5) == 2 + 0.5 * sqrt(2)
    end
    @testset "output never negative" begin
        @test RenewalDiD._approxcases(2, -1.5) == 0
    end
end

@testset "calculate numbers of infections in seed period" begin
    @test calcseedinfections1 == zeros(2, 2)
    @test calcseedinfections2 == zeros(3, 2)
    @test calcseedinfections3 == seedinfections1
    @testset for i in 1:2, j in 1:2 
        i == 2 && j == 1 && continue  # test for this value assumed error in Poisson approximation 
        @test calcseedinfections4[i, j] == predictedcalcseedinfections4[i, j]
    end
end 
        
@testset "calculate numbers of infections in study period" begin
    @test calcinfections1 == zeros(5, 2)
    @test calcinfections2 == zeros(6, 2)
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 3)), zeros(2, 2), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 3), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(3), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(4, 2)), zeros(2, 2), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(2), 3
    ) 
    @test calcinfections3 == predictedinfections3
    @test calcinfections4 == predictedinfections4
    @test calcinfections5[1, :] == predictedinfections5[1, :]  # test assumed error in Poisson approximation  
end     

@testset "proportion susceptible" begin
    infn = [0+1im  20+0.98im; 0+1im  40+0.96im; 20+0.8im  100+0.86im]
    expectedoutput = [1  1; 1  0.98; 1  0.96; 0.8  0.86]
    @testset for j in 1:2, t in 1:4
        @test RenewalDiD._prevpropsus(infn, t, j) == expectedoutput[t, j]
        @test (@inferred RenewalDiD._prevpropsus(infn, t, j)) == 
            RenewalDiD._prevpropsus(infn, t, j)    
    end
    @test_throws BoundsError RenewalDiD._prevpropsus(infn, 0, 1)
    @test_throws BoundsError RenewalDiD._prevpropsus(infn, 5, 2)
    @test_throws BoundsError RenewalDiD._prevpropsus(infn, 2, 3)
end

@testset "expected number of infections with proportion susceptible" begin
    # using `log(0)` and `log(1)` to emphasize that these parameters are natural logarithms
    @test RenewalDiD._expectedinfections(g_covid, 1, log(0), zeros(10)) == 0
    @test RenewalDiD._expectedinfections(g_covid, 1, log(1), ones(10)) == 
        sum(RenewalDiD.COVIDSERIALINTERVAL[2:11])
    @test RenewalDiD._expectedinfections(g_covid, 0, log(1), ones(10)) == 0
    @test RenewalDiD._expectedinfections(g_covid, 0.5, log(1), ones(10)) == 
        sum(RenewalDiD.COVIDSERIALINTERVAL[2:11]) / 2
    
    # `RenewalDiD._expectedinfections` now accepts proportions <0 and >1
end   

@testset "calculate numbers of infections in seed period with proportion susceptible" begin
    @test calcseedinfections5 == ones(2, 2) * (0 + 1im)
    @test calcseedinfections6 == ones(2, 2) * (0 + 1im)
    @test calcseedinfections7 == ones(3, 2) * (0 + 1im)
    @test calcseedinfections8 == predictedcalcseedinfections8
    @test calcseedinfections9 == predictedcalcseedinfections9
    @test calcseedinfections10[1, :] == predictedcalcseedinfections10[1, :]  # test assumed error in Poisson approximation
end 

@testset "calculate numbers of infections in study period with proportion susceptible" begin
    @test calcinfections6 == ones(5, 2) * (0 + 1im)
    @test calcinfections7 == ones(6, 2) * (0 + 1im)
        @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 3)), zeros(2, 2), ones(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 3), ones(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), ones(3), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(4, 2)), zeros(2, 2), ones(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, ComplexF64, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), ones(2), 3
    ) 
    @test calcinfections8 == predictedinfections8
    @testset for j in 1:2 
        @testset for t in 1:5 
            @test calcinfections9[t, j] ≈ predictedinfections9[t, j] atol=1e-9
        end
    end
end     

@testset "when tracking susceptibles, limit possible number of infections (seed)" begin
    infn = let 
        infn = zeros(ComplexF64, 6, 3)
        Mx = (m = zeros(6, 3); m[:, 2] .= -1.5; m[:, 3] .= 5; m)
        exptdseedcases = 10 .* ones(5, 3)
        Ns = 100 * ones(Int, 3)
        RenewalDiD._infections_seed!(infn, Mx, exptdseedcases, Ns, 5) 
        infn
    end
    @testset for t in 1:5 
        @test infn[t, 1] == 10 + (1 - 0.1 * t) * im
        # previous test assumed error in Poisson approximation
    end
    @testset for g in 1:3 
        @test minimum(real.(infn[:, g])) >= 0
        @test maximum(real.(infn[:, g])) <= 100
        @test sum(real.(infn[:, g])) <= 100
        @test minimum(imag.(infn[:, g])) >= 0
        @test maximum(imag.(infn[:, g])) <= 1
    end
end

@testset "when tracking susceptibles, limit possible infections (transmissions)" begin
    infn = let 
        _g(t) = generationtime(t; vec=(0.2 * ones(5)))
        infn = (m = zeros(ComplexF64, 6, 3); m[1, 1] = m[1, 2] = m[1, 3] = 1+0.99im; m)
        Mx = (m = zeros(6, 3); m[:, 2] .= -1.5; m[:, 3] .= 5; m)
        logR_0 = log(5) .* ones(5, 3)
        Ns = 100 * ones(Int, 3)
        RenewalDiD._infections_transmitted!(_g, infn, Mx, logR_0, Ns, 1)
        infn
    end
    @testset for g in 1:3 
        @test infn[1, g] == 1+0.99im
    end 
    @test infn[2, 1] ≈ 0.99 + (0.99^2) * im atol=1e-9
    @testset for t in 2:5 
        @test infn[t, 2] == 0+0.99im 
    end
    @testset for g in 1:3 
        @test minimum(real.(infn[:, g])) >= 0
        @test maximum(real.(infn[:, g])) <= 100
        @test sum(real.(infn[:, g])) <= 100
        @test minimum(imag.(infn[:, g])) >= 0
        @test maximum(imag.(infn[:, g])) <= 0.99
    end
end

@testset "output from `Base.show`" begin
    @test repr("text/plain", zerosstruct) == nzoutput 
    @test repr("text/plain", namedzerosstruct) == nnzoutput 
    @test repr("text/plain", namedzerostructunlimitied) == unnzoutput 
end

@testset "expected variance of _approxcasescalc" begin
    rng = Xoshiro(1)
    x = 50 
    sigmas = rand(rng, Normal(0, 1), 100_000)
    vecofvalues = [RenewalDiD._approxcasescalc(x, sigma) for sigma in sigmas] 
    @test mean(vecofvalues) ≈ 50 rtol=1e-2
    @test var(vecofvalues) ≈ 50 rtol=1e-2
end
