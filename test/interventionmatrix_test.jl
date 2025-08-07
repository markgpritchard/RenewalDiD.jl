# test functions related to generating and using `InterventionMatrix` structs

using RenewalDiD
using Test

M1 = InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true) 
M2 = InterventionMatrix{Int}(3, [2, nothing])
M3 = InterventionMatrix(10, 4, 5, 8, 8, nothing)
M4 = InterventionMatrix{Float64}(5, 2, 3, nothing)   

M1expectedcontents = [
    0  0
    1  0 
    1  1
]
M1expectedtimestring = ["1", "2", "3"]
M1expectedcontentsstring = [
    "0"  "0"
    "1"  "0" 
    "1"  "1"
]
M1expectedcombinedstring = [
    "1"  "0"  "0"
    "2"  "1"  "0" 
    "3"  "1"  "1"
]
M1expectedheaderstring = ["time", "1", "2"]
M1expectedshowoutput = "3×2 InterventionMatrix{Int64}\n time │ 1  2\n──────┼──────\n    1 \
    │ 0  0\n    2 │ 1  0\n    3 │ 1  1\n"
M3expectedcontents = [
    0  0  0  0  0
    1  0  0  0  0
    1  1  0  0  0
    1  1  1  1  0
    1  1  1  1  0
]
M3expectedtimestring = ["1", "⋮", "4", "5", "⋮", "8", "⋮", "10"]
M3expectedcontentsstring = [
    "0"  "0"  "0"  "0"  "0"
    "⋮"  "⋮"  "⋮"  "⋮"  "⋮"
    "1"  "0"  "0"  "0"  "0"
    "1"  "1"  "0"  "0"  "0"
    "⋮"  "⋮"  "⋮"  "⋮"  "⋮"
    "1"  "1"  "1"  "1"  "0"
    "⋮"  "⋮"  "⋮"  "⋮"  "⋮"
    "1"  "1"  "1"  "1"  "0"
]
M4expectedcontentsstring = [
    "0.0"  "0.0"  "0.0" 
    "1.0"  "0.0"  "0.0" 
    "1.0"  "1.0"  "0.0"
    "⋮"  "⋮"  "⋮"
    "1.0"  "1.0"  "0.0"
]
M4expectedcombinedstring = [
    "1"  "0.0"  "0.0"  "0.0" 
    "2"  "1.0"  "0.0"  "0.0" 
    "3"  "1.0"  "1.0"  "0.0"
    "⋮"  "⋮"  "⋮"  "⋮"
    "5"  "1.0"  "1.0"  "0.0"
]
M4expectedheaderstring = ["time", "1", "2", "3"]
M4expectedshowoutput ="5×3 InterventionMatrix{Float64}\n time │   1    2    3\n──────┼─────\
    ──────────\n    1 │ 0.0  0.0  0.0\n    2 │ 1.0  0.0  0.0\n    3 │ 1.0  1.0  0.0\n    ⋮ \
    │   ⋮    ⋮    ⋮\n    5 │ 1.0  1.0  0.0\n"
tw1 = "InterventionMatrix with no intervention in any group"
tw2 = "All groups in InterventionMatrix have intervention before end of duration"

@testset "constructing an InterventionMatrix" begin
    @test InterventionMatrix{Int}(2, 2; mutewarnings=true) isa AbstractMatrix
    # behaviour of `mutewarnings` is tested below
    @test size(InterventionMatrix{Int}(2, 2; mutewarnings=true)) == (2, 1)
    @test size(InterventionMatrix{Int}(3, 2; mutewarnings=true)) == (3, 1)
    @test size(InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true)) == (3, 2)
end
@testset "construction without an explicit type" begin
    # test that contructor without a type works, not that it  generates a type `Int`; tests 
    # of `Base.show` below demonstrate that the default type is `Int`
    _M1a = InterventionMatrix(2, 2; mutewarnings=true)
    _M1b = InterventionMatrix{Int}(2, 2; mutewarnings=true) 
    _M2a = InterventionMatrix(3, 2; mutewarnings=true)
    _M2b = InterventionMatrix{Int}(3, 2; mutewarnings=true)
    @test _M1a == _M1b
    @test _M2a == _M2b
end
@testset "construction with different arrangement of arguments" begin
    _M1a = InterventionMatrix(4, 2, 3; mutewarnings=true)
    _M1b = InterventionMatrix{Int}(4, [2, 3]; mutewarnings=true)
    _M2a = InterventionMatrix(4, 2, 3, 4; mutewarnings=true)
    _M2b = InterventionMatrix{Int}(4, [2, 3, 4]; mutewarnings=true)
    @test InterventionMatrix(3, [2, 3]; mutewarnings=true) == M1
    @test InterventionMatrix(3, 2, 3; mutewarnings=true) == M1
    @test _M1a == _M1b
    @test _M2a == _M2b
end
@testset "`starttimes` greater that `duration` equivalent to `nothing`" begin
    @test InterventionMatrix(4, 2, 3, nothing) == InterventionMatrix{Int}(4, [2, 3, 5])
end
@testset "warnings during construction" begin
    @test_warn tw1 InterventionMatrix(2, nothing, nothing)
    @test_nowarn InterventionMatrix{Int}(2, [2, 3])
    @test_nowarn InterventionMatrix(4, 2, 3, nothing)
    @test_nowarn InterventionMatrix(2, nothing, nothing; mutewarnings=true)
    @test_warn tw1 InterventionMatrix(2, nothing, nothing; mutewarnings=false)
    @test_warn tw2 InterventionMatrix(3, 2, 3)
    @test_nowarn InterventionMatrix(3, 2, 3; mutewarnings=true)
    @test_warn tw2 InterventionMatrix(3, 2, 3; mutewarnings=false)
end
@testset "construction with explicit non-`Int` type" begin
    @test InterventionMatrix{Bool}(3, [2, 3]; mutewarnings=true)[3, 1]  
end
@testset "invalid arguments to contructor" begin
    @test_throws InexactError InterventionMatrix(3.3, 2; mutewarnings=true)  # a change in 
    # the function means that it assesses the timings of interventions before finding the 
    # InexactError
    @test_throws InexactError InterventionMatrix(3, 2.3)
    @test_throws InexactError InterventionMatrix(3, [2.3, 2])
end
@testset "indexing InterventionMatrix" begin
    @test M1[1] == 0
    @test M1[1, 1] == 0
    @test M1[3, 1] == 1
    @test M1[3] == 1
    @test M1[4] == 0
    @test M1[6] == 1
    @test_throws BoundsError M1[7]
    @test_throws BoundsError M1[4, 1]
    @test_throws BoundsError M1[1, 3]
end
@testset "indexing InterventionMatrix when `starttimes` are missing" begin
    @test M2[1, 2] == 0
    @test M2[3, 2] == 0
end
@testset "`_duration`" begin
    @test RenewalDiD._duration(M1) == 3
    @test RenewalDiD._duration(InterventionMatrix(4, [2, 3]; mutewarnings=true)) == 4
end
@testset "`_showtimes`" begin
    @test RenewalDiD._showtimes(M1) == [1, 2, 3]
    @test RenewalDiD._showtimes(M3) == [1, 4, 5, 8, 10]
end
@testset "`_showcontents`" begin
    @test RenewalDiD._showcontents(M1, RenewalDiD._showtimes(M1)) == M1expectedcontents
    @test RenewalDiD._showcontents(M3, RenewalDiD._showtimes(M3)) == M3expectedcontents
end
@testset "`_showstrings`" begin
    ts1, cs1 = RenewalDiD._showstrings(M1)
    ts3, cs3 = RenewalDiD._showstrings(M3)
    ts4, cs4 = RenewalDiD._showstrings(M4)
    @test ts1 == M1expectedtimestring
    @test cs1 == M1expectedcontentsstring 
    @test ts3 == M3expectedtimestring
    @test cs3 == M3expectedcontentsstring
    @test cs4 == M4expectedcontentsstring
end
@testset "`_showcombinedstrings`" begin
    @test RenewalDiD._showcombinedstrings(M1) == M1expectedcombinedstring
    @test RenewalDiD._showcombinedstrings(M4) == M4expectedcombinedstring
end
@testset "`_showheader`" begin
    @test RenewalDiD._showheader(M1) == M1expectedheaderstring
    @test RenewalDiD._showheader(M4) == M4expectedheaderstring
end
@testset "output from `Base.show`" begin
    @test repr("text/plain", M1) == M1expectedshowoutput
    @test repr("text/plain", M4) == M4expectedshowoutput
end
