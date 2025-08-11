# test functions related to generating and using offset intervention structs

using RenewalDiD
using Test

IM1 = InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true) 
IM2 = InterventionMatrix{Int}(3, [2, nothing])
IM3 = InterventionMatrix(10, 4, 5, 8, 8, nothing)
IM4 = InterventionMatrix{Float64}(5, 2, 3, nothing)   
OIM1 = OffsetInterventionMatrix{Int}(3, 0, [2, 3]; mutewarnings=true) 
OIM2 = OffsetInterventionMatrix{Int}(3, 0, [2, nothing])
OIM3 = OffsetInterventionMatrix(10, 0, 4, 5, 8, 8, nothing)
OIM4 = OffsetInterventionMatrix{Float64}(5, 0, 2, 3, nothing)  
OIM5 = OffsetInterventionMatrix{Int}(3, -5, [7, 8]; mutewarnings=true) 

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
M1expectedshowoutput = "3×2 OffsetInterventionMatrix{Int64} (offset 0)\n time │ 1  \
    2\n──────┼──────\n    1 │ 0  0\n    2 │ 1  0\n    3 │ 1  1\n"
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
M4expectedshowoutput ="5×3 OffsetInterventionMatrix{Float64} (offset 0)\n time │   1    \
    2    3\n──────┼───────────────\n    1 │ 0.0  0.0  0.0\n    2 │ 1.0  0.0  0.0\n    3 │ \
    1.0  1.0  0.0\n    ⋮ │   ⋮    ⋮    ⋮\n    5 │ 1.0  1.0  0.0\n"
M5expectedshowoutput ="3×2 OffsetInterventionMatrix{Int64} (offset -5)\n time │ 1  \
    2\n──────┼──────\n    1 │ 0  0\n    2 │ 1  0\n    3 │ 1  1\n"
tw1 = "OffsetInterventionMatrix with no intervention in any group"
tw2 = "All groups in OffsetInterventionMatrix have intervention before end of duration"

@testset "constructing an OffsetInterventionMatrix" begin
    @test OffsetInterventionMatrix{Int}(2, 0, 2; mutewarnings=true) isa AbstractMatrix
    # behaviour of `mutewarnings` is tested below
    @test size(OffsetInterventionMatrix{Int}(2, 0, 2; mutewarnings=true)) == (2, 1)
    @test size(OffsetInterventionMatrix{Int}(3, 0, 2; mutewarnings=true)) == (3, 1)
    @test size(OffsetInterventionMatrix{Int}(3, 0, [2, 3]; mutewarnings=true)) == (3, 2)
    @test size(OffsetInterventionMatrix{Int}(2, 1, 2; mutewarnings=true)) == (2, 1)
    @test size(OffsetInterventionMatrix{Int}(3, 1, 2; mutewarnings=true)) == (3, 1)
    @test size(OffsetInterventionMatrix{Int}(3, 1, [2, 3]; mutewarnings=true)) == (3, 2)
    @test size(OffsetInterventionMatrix{Int}(2, -1, 2; mutewarnings=true)) == (2, 1)
    @test size(OffsetInterventionMatrix{Int}(3, -1, 2; mutewarnings=true)) == (3, 1)
    @test size(OffsetInterventionMatrix{Int}(3, -1, [2, 3]; mutewarnings=true)) == (3, 2)
end
@testset "construction without an explicit type" begin
    # test that contructor without a type works, not that it  generates a type `Int`; tests 
    # of `Base.show` below demonstrate that the default type is `Int`
    _M1a = OffsetInterventionMatrix(2, -1, 2; mutewarnings=true)
    _M1b = OffsetInterventionMatrix{Int}(2, -1, 2; mutewarnings=true) 
    _M2a = OffsetInterventionMatrix(3, -1, 2; mutewarnings=true)
    _M2b = OffsetInterventionMatrix{Int}(3, -1, 2; mutewarnings=true)
    @test _M1a == _M1b
    @test _M2a == _M2b
end
@testset "construction with different arrangement of arguments" begin
    _M1a = OffsetInterventionMatrix(4, -1, 2, 3; mutewarnings=true)
    _M1b = OffsetInterventionMatrix{Int}(4, -1, [2, 3]; mutewarnings=true)
    _M2a = OffsetInterventionMatrix(4, -1, 2, 3, 4; mutewarnings=true)
    _M2b = OffsetInterventionMatrix{Int}(4, -1, [2, 3, 4]; mutewarnings=true)
    @test OffsetInterventionMatrix(3, 0, [2, 3]; mutewarnings=true) == OIM1
    @test OffsetInterventionMatrix(3, 0, 2, 3; mutewarnings=true) == OIM1
    @test _M1a == _M1b
    @test _M2a == _M2b
end
@testset "offset of zero equal to InterventionMatrix" begin
    @test IM1 == OIM1 
    @test IM2 == OIM2 
    @test IM3 == OIM3 
    @test IM4 == OIM4 
end
@testset "`starttimes` greater that `duration` equivalent to `nothing`" begin
    @test OffsetInterventionMatrix(5, 0, 2, 3, nothing) == OIM4
    @test OffsetInterventionMatrix(5, 0, 2, 3, 6) == OIM4
    @test OffsetInterventionMatrix(5, 0, 2, 3, 0) == OIM4
    @test OffsetInterventionMatrix(5, 1, 1, 2, nothing) == OIM4
    @test OffsetInterventionMatrix(5, 1, 1, 2, 5) == OIM4
    @test OffsetInterventionMatrix(5, 1, 1, 2, 0) == OIM4
    @test OffsetInterventionMatrix(5, 1, 1, 2, 1; mutewarnings=true) != OIM4
end
@testset "retain raw data in subsequent matrices" begin
    _IM = InterventionMatrix(10, 1, 5, 9, 11, 15, nothing)
    _OM1 = OffsetInterventionMatrix(_IM, 0)
    _OM2 = OffsetInterventionMatrix(_OM1, 0)
    _OM3 = OffsetInterventionMatrix(_OM2, -1)
    _OM4 = OffsetInterventionMatrix(_OM3, 1)
    _OM5 = OffsetInterventionMatrix(_OM4, 1)
    _OM6 = OffsetInterventionMatrix(_OM5, -9)
    @test _IM == _OM1
    @test _OM1 == _OM2
    @test _OM3 == InterventionMatrix(10, nothing, 4, 8, 10, nothing, nothing)
    @test _OM4 == _IM
    @test _OM5 == InterventionMatrix(10, 2, 6, 10, nothing, nothing, nothing)
    @test _OM6 == InterventionMatrix(10, nothing, nothing, 1, 3, 7, nothing)
end
@testset "warnings and errors during construction" begin
    @test_throws MethodError OffsetInterventionMatrix(2, nothing, nothing)
    @test_throws MethodError OffsetInterventionMatrix(2, nothing, [2, 3])
    @test_throws MethodError OffsetInterventionMatrix(2, [2, 3])
    @test_warn tw1 OffsetInterventionMatrix(2, 0, nothing, nothing)
    @test_nowarn OffsetInterventionMatrix{Int}(2, 0, [2, 3])
    @test_nowarn OffsetInterventionMatrix(4, 0, 2, 3, nothing)
    @test_nowarn OffsetInterventionMatrix(2, 0, nothing, nothing; mutewarnings=true)
    @test_warn tw1 OffsetInterventionMatrix(2, 0, nothing, nothing; mutewarnings=false)
    @test_warn tw2 OffsetInterventionMatrix(3, 0, 2, 3)
    @test_nowarn OffsetInterventionMatrix(3, 0, 2, 3; mutewarnings=true)
    @test_warn tw2 OffsetInterventionMatrix(3, 0, 2, 3; mutewarnings=false)
end
@testset "construction with explicit non-`Int` type" begin
    @test OffsetInterventionMatrix{Bool}(3, 0, [2, 4])[3, 1]  
    @test !OffsetInterventionMatrix{Bool}(3, 0, [2, 4])[3, 2]  
    @test OffsetInterventionMatrix{Bool}(3, -1, [2, 4])[3, 2]  
end
@testset "invalid arguments to contructor" begin
    @test_throws InexactError OffsetInterventionMatrix(3.3, 0, 2; mutewarnings=true)
    @test_throws InexactError OffsetInterventionMatrix(3, 0.5, 2)
    @test_throws InexactError OffsetInterventionMatrix(3, 0, 2.3)
    @test_throws InexactError OffsetInterventionMatrix(3, 0, [2.3, 2])
    @test_throws ArgumentError OffsetInterventionMatrix(3, 0, 1, [2, 3])
    @test_throws ArgumentError OffsetInterventionMatrix(3, 0, [2, 3], 4)
end
@testset "indexing InterventionMatrix" begin
    @testset for M in [OIM1, OIM5]
        @test M[1] == 0
        @test M[1, 1] == 0
        @test M[3, 1] == 1
        @test M[3] == 1
        @test M[4] == 0
        @test M[6] == 1
        @test_throws BoundsError M[7]
        @test_throws BoundsError M[4, 1]
        @test_throws BoundsError M[1, 3]
    end
end
@testset "indexing InterventionMatrix when `starttimes` are missing" begin
    @test OIM2[1, 2] == 0
    @test OIM2[3, 2] == 0
end
@testset "`_duration`" begin
    @test RenewalDiD._duration(OIM1) == 3
    @test RenewalDiD._duration(OffsetInterventionMatrix(4, 2, [2, 3]; mutewarnings=true)) == 4
end
@testset "`_showtimes`" begin
    @test RenewalDiD._showtimes(OIM1) == [1, 2, 3]
    @test RenewalDiD._showtimes(OIM3) == [1, 4, 5, 8, 10]
    @test RenewalDiD._showtimes(OffsetInterventionMatrix(OIM3, 5)) == [1, 9, 10]
end
@testset "`_showcontents`" begin
    @test RenewalDiD._showcontents(OIM1, RenewalDiD._showtimes(OIM1)) == M1expectedcontents
    @test RenewalDiD._showcontents(OIM3, RenewalDiD._showtimes(OIM3)) == M3expectedcontents
    @test RenewalDiD._showcontents(OIM5, RenewalDiD._showtimes(OIM5)) == M1expectedcontents
end
@testset "`_showstrings`" begin
    ts1, cs1 = RenewalDiD._showstrings(OIM1)
    ts2, cs2 = RenewalDiD._showstrings(OIM5)
    ts3, cs3 = RenewalDiD._showstrings(OIM3)
    ts4, cs4 = RenewalDiD._showstrings(OIM4)
    @test ts1 == M1expectedtimestring
    @test cs1 == M1expectedcontentsstring     
    @test ts2 == M1expectedtimestring
    @test cs2 == M1expectedcontentsstring 
    @test ts3 == M3expectedtimestring
    @test cs3 == M3expectedcontentsstring
    @test cs4 == M4expectedcontentsstring
end
@testset "`_showcombinedstrings`" begin
    @test RenewalDiD._showcombinedstrings(OIM1) == M1expectedcombinedstring
    @test RenewalDiD._showcombinedstrings(OIM5) == M1expectedcombinedstring
    @test RenewalDiD._showcombinedstrings(OIM4) == M4expectedcombinedstring
end
@testset "`_showheader`" begin
    @test RenewalDiD._showheader(OIM1) == M1expectedheaderstring
    @test RenewalDiD._showheader(OIM5) == M1expectedheaderstring
    @test RenewalDiD._showheader(OIM4) == M4expectedheaderstring
end
@testset "output from `Base.show`" begin
    @test repr("text/plain", OIM1) == M1expectedshowoutput
    @test repr("text/plain", OIM4) == M4expectedshowoutput
    @test repr("text/plain", OIM5) == M5expectedshowoutput
end
