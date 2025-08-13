# test functions related to generating and using `InterventionMatrix` structs

using RenewalDiD
using Test

M1 = InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true) 
M2 = InterventionMatrix{Int}(3, [2, nothing])
M3 = InterventionMatrix(10, 4, 5, 8, 8, nothing)
M4 = InterventionMatrix{Float64}(5, 2, 3, nothing)   
M5 = InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true, offset=0) 
M6 = InterventionMatrix{Int}(3, [2, nothing]; offset=0)
M7 = InterventionMatrix(10, 4, 5, 8, 8, nothing; offset=0)
M8 = InterventionMatrix{Float64}(5, 2, 3, nothing; offset=0)  
M9 = InterventionMatrix{Int}(3, [7, 8]; mutewarnings=true, offset=-5) 

M1expectedcontents = [0  0; 1  0; 1  1]
M1expectedtimestring = ["1", "2", "3"]
M1expectedcontentsstring = ["0"  "0"; "1"  "0"; "1"  "1"]
M1expectedcombinedstring = ["1"  "0"  "0"; "2"  "1"  "0"; "3"  "1"  "1"]
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
    @test InterventionMatrix{Int}(2, 2; mutewarnings=true) isa AbstractInterventionArray
    # (behaviour of `mutewarnings` is tested below)
    @test InterventionMatrix{Int}(2, 2; mutewarnings=true) isa AbstractMatrix
    @test size(InterventionMatrix{Int}(2, 2; mutewarnings=true)) == (2, 1)
    @test size(InterventionMatrix{Int}(3, 2; mutewarnings=true)) == (3, 1)
    @test size(InterventionMatrix{Int}(3, [2, 3]; mutewarnings=true)) == (3, 2)
end
@testset "construction without an explicit type" begin
    # test that contructor without a type works, not that it  generates a type `Int`; later
    # tests of `Base.show` demonstrate that the default type is `Int`
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
@testset "offset of zero equal to InterventionMatrix" begin
    @test M1 == M5 
    @test M2 == M6 
    @test M3 == M7 
    @test M4 == M8 
end
@testset "`starttimes` greater that `duration` or equal to 1 equivalent to `nothing`" begin
    @test InterventionMatrix{Int}(5, 2, 3, nothing) == M4
    @test InterventionMatrix{Int}(5, 2, 3, 6) == M4
    @test InterventionMatrix{Int}(5, 2, 3, 1) == M4
    @test InterventionMatrix{Int}(5, 2, 3, -10) == M4
    @test InterventionMatrix{Int}(5, 1, 2, nothing; offset=1) == M4
    @test InterventionMatrix{Int}(5, 1, 2, 5; offset=1) == M4
    @test InterventionMatrix{Int}(5, 1, 2, 0; offset=1) == M4
    @test InterventionMatrix{Int}(5, 1, 2, 1; mutewarnings=true, offset=1) != M4
end
@testset "retain offset and raw data in subsequent matrices" begin
    _IM = InterventionMatrix(10, 1, 5, 9, 11, 15, nothing)
    _OM1 = InterventionMatrix(_IM; offset=0)
    _OM2 = InterventionMatrix(_OM1; offset=0)
    _OM3 = InterventionMatrix(_OM2; offset=-1)
    _OM4 = InterventionMatrix(_OM3; offset=1)
    _OM5 = InterventionMatrix(_OM4; offset=1)
    _OM6 = InterventionMatrix(_OM5; offset=-9)
    @test _IM == _OM1
    @test _OM1 == _OM2
    @test _OM3 == InterventionMatrix(10, nothing, 4, 8, 10, nothing, nothing)
    @test _OM4 == _IM
    @test _OM5 == InterventionMatrix(10, 2, 6, 10, nothing, nothing, nothing)
    @test _OM6 == InterventionMatrix(10, nothing, nothing, 1, 3, 7, nothing)
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
    @test InterventionMatrix{Bool}(3, [2, 4])[3, 1]  
    @test !InterventionMatrix{Bool}(3, [2, 4])[3, 2]  
    @test InterventionMatrix{Bool}(3, [2, 4]; offset=-1)[3, 2]  
end
@testset "invalid arguments to contructor" begin
    @test_throws InexactError InterventionMatrix(3.3, 2; mutewarnings=true)  # a change in 
    # the function means that it assesses the timings of interventions before finding the 
    # InexactError
    @test_throws InexactError InterventionMatrix(3, 2.3)
    @test_throws InexactError InterventionMatrix(3, [2.3, 2])
    @test_throws InexactError InterventionMatrix(3, 2; offset=0.5)
    @test_throws InexactError InterventionMatrix(3, [2, 2]; offset=0.5)
    @test_throws ArgumentError InterventionMatrix(3, 1, [2, 3])
    @test_throws ArgumentError InterventionMatrix(3, [2, 3], 4)
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
    @test RenewalDiD._showtimes(M1, nothing) == [1, 2, 3]
    @test RenewalDiD._showtimes(M3, nothing) == [1, 4, 5, 8, 10]
end
@testset "`_showcontents`" begin
    @test RenewalDiD._showcontents(M1, RenewalDiD._showtimes(M1, nothing), nothing) == M1expectedcontents
    @test RenewalDiD._showcontents(M3, RenewalDiD._showtimes(M3, nothing), nothing) == M3expectedcontents
end
@testset "`_showstrings`" begin
    ts1, cs1 = RenewalDiD._showstrings(M1, nothing)
    ts3, cs3 = RenewalDiD._showstrings(M3, nothing)
    ts4, cs4 = RenewalDiD._showstrings(M4, nothing)
    @test ts1 == M1expectedtimestring
    @test cs1 == M1expectedcontentsstring 
    @test ts3 == M3expectedtimestring
    @test cs3 == M3expectedcontentsstring
    @test cs4 == M4expectedcontentsstring
end
@testset "`_showcombinedstrings`" begin
    @test RenewalDiD._showcombinedstrings(M1, nothing) == M1expectedcombinedstring
    @test RenewalDiD._showcombinedstrings(M4, nothing) == M4expectedcombinedstring
end
@testset "`_showheader`" begin
    @test RenewalDiD._showheader(M1) == M1expectedheaderstring
    @test RenewalDiD._showheader(M4) == M4expectedheaderstring
end
@testset "output from `Base.show`" begin
    @test repr("text/plain", M1) == M1expectedshowoutput
    @test repr("text/plain", M4) == M4expectedshowoutput
end
@testset "concatenate vectors and matrices" begin
    v1 = InterventionVector(5, 2)
    v2 = InterventionVector(5, nothing)
    v3 = InterventionVector(6, nothing)
    _M1 = InterventionMatrix(5, 2; mutewarnings=true)
    _M2 = InterventionMatrix(5, nothing; mutewarnings=true)
    _M3 = InterventionMatrix(6, [3, nothing])
    _M4 = InterventionMatrix(5, [3, nothing])
    M1 = cat(v1; dims=2, mutewarnings=true)
    M1a = cat(v1; dims=Val{2}(), mutewarnings=true)
    M2 = cat(v1, v2; dims=2)
    M3 = cat(M1; dims=2, mutewarnings=true)
    M4 = cat(M2; dims=2)
    M5 = cat(_M1, v2; dims=2)
    M5a = cat(v1, _M2; dims=2)
    M5b = cat(_M1, _M2; dims=2)
    M6 = cat(v3, _M3; dims=2)
    M7 = cat(_M4, M2; dims=2)

    @test_throws ArgumentError cat(v1; dims=1)
    @test_throws ArgumentError cat(v1; dims=Val{1}())
    @test_throws ArgumentError cat(v1, v2; dims=1)
    @test_throws ArgumentError cat(M1; dims=1)
    @test_throws ArgumentError cat(M2; dims=1)
    @test_throws ArgumentError cat(_M1, v2; dims=1)
    @test_throws ArgumentError cat(v1, _M2; dims=1)
    @test_throws ArgumentError cat(_M1, _M2; dims=1)
    @test_throws ArgumentError cat(v3, _M3; dims=1)
    @test_throws ArgumentError cat(_M4, M2; dims=1)

    @test_warn tw2 cat(v1; dims=2)
    @test_nowarn cat(v1; dims=2, mutewarnings=true)
    @test_warn tw2 cat(v1; dims=2, mutewarnings=false)
    @test_nowarn cat(v1, v2; dims=2)
    @test_warn tw2 cat(M1; dims=2)
    @test_nowarn cat(M1; dims=2, mutewarnings=true)
    @test_warn tw2 cat(M1; dims=2, mutewarnings=false)
    @test_nowarn cat(M2; dims=2)

    @test M1 == _M1
    @test M1a == _M1
    @test v1 != _M1
    @test M2 == InterventionMatrix(5, [2, nothing])
    @test_throws DimensionMismatch cat(v1, v3; dims=2)
    @test M3 == M1
    @test M4 == M2
    @test M5 == M2
    @test M5a == M2
    @test M5b == M2
    @test M6 == InterventionMatrix(6, [nothing, 3, nothing])
    @test_throws DimensionMismatch cat(v1, _M3; dims=2)
    @test_throws DimensionMismatch cat(M1, _M3; dims=2)
    @test_throws DimensionMismatch cat(M2, _M3; dims=2)
    @test M7 == InterventionMatrix(5, [3, nothing, 2, nothing])
end
