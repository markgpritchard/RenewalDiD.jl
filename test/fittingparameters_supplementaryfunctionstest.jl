# test functions that aid `RenewalDiD._renewaldid`

using RenewalDiD
using Test

M1 = InterventionMatrix(4, [2, 3, nothing]) 
M2 = InterventionMatrix(5, [nothing, nothing]; mutewarnings=true) 

obs1 = (_obs = zeros(Int, 20, 2); _obs[1, :] .+= 1; _obs[1:5, :] .+= 1; _obs)
obs2 = (_obs = zeros(Int, 20, 2); _obs[1, 1] += 1; _obs[1:5, :] .+= 1; _obs)

seedinfections1 = [1  0; 2  1]
M_x1 = [2  0; -0.5  1; 1  0.5; -2  1; 0  1]
calcinfections1 = RenewalDiD._infections(
    g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(2), 2
)
calcinfections2 = RenewalDiD._infections(
    g_covid, zeros(6, 2), log.(zeros(4, 2)), zeros(2, 2), zeros(2), 2
)
calcinfections3 = RenewalDiD._infections(
    generationtime, zeros(5, 2), log.(zeros(3, 2)), seedinfections1, 1000 .* ones(2), 2;
    vec=[0, 1],
)
calcinfections4 = RenewalDiD._infections(
    generationtime, zeros(5, 2), log.(1.5 .* ones(3, 2)), seedinfections1, 1000 .* ones(2), 2;
    vec=[0, 1],
)
calcinfections5 = RenewalDiD._infections(
    generationtime, M_x1, log.(ones(3, 2)), seedinfections1, 1000 .* ones(2), 2;
    vec=[0, 1]
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
predictedinfections3 = [1  0; 2  1; 0  0; 0  0; 0  0]
predictedinfections4 = [1  0; 2  1; 3  1.5; 4.5  2.25; 6.75  3.375]
predictedinfections5 = [3  0; 1  2; 2  3; 0  6; 0  12]

@testset "calculate ntimes" begin
    @test RenewalDiD._ntimes(zeros(2, 3)) == 2
    @test RenewalDiD._ntimes(zeros(3, 2)) == 3
    @test RenewalDiD._ntimes(M1) == 4
    @test_throws MethodError RenewalDiD._ntimes(zeros(2))
end

@testset "calculate ngroups" begin
    @test RenewalDiD._ngroups(zeros(2, 3)) == 3
    @test RenewalDiD._ngroups(zeros(3, 2)) == 2
    @test RenewalDiD._ngroups(M1) == 3
    @test_throws MethodError RenewalDiD._ngroups(zeros(2))
end

@testset "generate vector of gamma values" begin
    @test RenewalDiD._gammavec(0, 1, zeros(3)) == zeros(3)
    @test RenewalDiD._gammavec(0, 1, zeros(4)) == zeros(4)
    @test RenewalDiD._gammavec(1, 1, zeros(3)) == ones(3)
    @test RenewalDiD._gammavec(1, 1, [1, -1, 0]) == [2, 0, 1]
    @test RenewalDiD._gammavec(1, 0.5, [1, -1, 0]) == [1.5, 0.5, 1]
end

@testset "generate vector of theta values" begin
    @test RenewalDiD._thetavec(0, zeros(3), 1) == zeros(4)
    @test RenewalDiD._thetavec(0, zeros(4), 1) == zeros(5)
    @test RenewalDiD._thetavec(1, zeros(3), 1) == ones(4)
    @test RenewalDiD._thetavec(1, ones(3), 1) == [1, 2, 3, 4]
    @test RenewalDiD._thetavec(1, ones(3), 0.5) == [1, 1.5, 2, 2.5]
    @test RenewalDiD._thetavec(1, [-2.5, 0.5, 0.5], 0.5) == [1, -0.25, 0, 0.25]
end

@testset "generate matrix of log R_0 values" begin
    @test RenewalDiD._predictedlogR_0(0, zeros(3), zeros(4), 0, M1) == zeros(4, 3)
    @test RenewalDiD._predictedlogR_0(0, zeros(2), zeros(5), 2, M2) == zeros(5, 2)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(2), zeros(4), 0, M1)
    @test_throws DimensionMismatch RenewalDiD._predictedlogR_0(0, zeros(3), zeros(3), 0, M1)
    @test RenewalDiD._predictedlogR_0(1, zeros(3), zeros(4), 0, M1) == ones(4, 3)
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], zeros(4), 0, M1) == R0prediction1
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], 0, M1) == R0prediction2
    @test RenewalDiD._predictedlogR_0(1, [-1, 0, 1], [-1, 0.5, 2.5, 0], -1, M1) == R0prediction3
end

@testset "expected seed cases" begin
    @test RenewalDiD._expectedseedcases(zeros(20, 2), 7) == zeros(7, 2)
    @test RenewalDiD._expectedseedcases(zeros(20, 2), 5) == zeros(5, 2)
    @test RenewalDiD._expectedseedcases(zeros(20, 3), 5) == zeros(5, 3)
    @test RenewalDiD._expectedseedcases(obs1, 5) == seedexpectation1
    @test RenewalDiD._expectedseedcases(obs1, 5; doubletime=5) == seedexpectation2
    @test RenewalDiD._expectedseedcases(obs1, 5; doubletime=10) == seedexpectation3
    @test RenewalDiD._expectedseedcases(obs1, 5; sampletime=5) == seedexpectation4
    @test RenewalDiD._expectedseedcases(obs1, 5; sampletime=3) == seedexpectation5
    @test RenewalDiD._expectedseedcases(obs2, 5) == seedexpectation6
    @test RenewalDiD._expectedseedcases(zeros(2, 3), 5) == zeros(5, 3)
    @test_throws BoundsError RenewalDiD._expectedseedcases(zeros(2, 3), 5; sampletime=5)
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
    @test RenewalDiD._expectedinfections(g_seir, log(0), zeros(10); gamma=0.5) == 0
    @test RenewalDiD._expectedinfections(g_seir, log(1), [1, 0, 0, 0, 0, 0]; gamma=0.5) == 
        0.07468060255179591
    @test_throws UndefKeywordError RenewalDiD._expectedinfections(g_seir, log(0), zeros(10))
    # we do not need to use `generationtime` with `g_seir` but test it as though `g_seir` 
    # were a user-generated function
    @test RenewalDiD._expectedinfections(
        generationtime, log(0), zeros(10); 
        func=g_seir, gamma=0.5,
    ) == 0
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), [1, 0, 0, 0, 0, 0]; 
        func=g_seir, gamma=0.5,
    ) == 0.07468060255179591
    @test RenewalDiD._expectedinfections(
        generationtime, log(1), [1, 0, 0, 0, 0, 0]; 
        func=g_seir, gamma=0.5, t_max=5,
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
        func=g_seir, gamma=0.5, vec=[0, 0, 0, 0, 0, 1],
    ) 
end   
        
@testset "approximate number of cases" begin
    @testset "with no variance is just the point estimate" begin
        @test RenewalDiD._approxcases(0, 0) == 0
        @test RenewalDiD._approxcases(1, 0) == 1
    end
    @testset "account for variance" begin
        @test RenewalDiD._approxcases(1, 0.5) == 1.5
        @test RenewalDiD._approxcases(2, 0.5) == 3
    end
    @testset "output never negative" begin
        @test RenewalDiD._approxcases(2, -1.5) == 0
    end
end
        
@testset "calculate numbers of infections" begin
    @test calcinfections1 == zeros(5, 2)
    @test calcinfections2 == zeros(6, 2)
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 3)), zeros(2, 2), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 3), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(3), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(4, 2)), zeros(2, 2), zeros(2), 2
    ) 
    @test_throws DimensionMismatch RenewalDiD._infections(
        g_covid, zeros(5, 2), log.(zeros(3, 2)), zeros(2, 2), zeros(2), 3
    ) 
    @test calcinfections3 == predictedinfections3
    @test calcinfections4 == predictedinfections4
    @test calcinfections5 == predictedinfections5
end     
