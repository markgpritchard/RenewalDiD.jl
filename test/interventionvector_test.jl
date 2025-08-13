# test functions related to generating and using `InterventionMatrix` structs

using RenewalDiD
using Test

v1 = InterventionVector{Int}(10, 5) 
v2 = InterventionVector{Int}(15, nothing) 

v1expectedshowoutput = "10-element InterventionVector{Int64}\n time │ 1\n──────┼───\n    1 \
    │ 0\n    ⋮ │ ⋮\n    5 │ 1\n    ⋮ │ ⋮\n   10 │ 1\n"
v2expectedshowoutput = "15-element InterventionVector{Int64}\n time │ 1\n──────┼───\n    1 \
    │ 0\n    ⋮ │ ⋮\n   15 │ 0\n"

@testset "constructing an InterventionMatrix" begin
    @test InterventionVector{Int}(2, 2) isa AbstractInterventionArray
    # (behaviour of `mutewarnings` is tested below)
    @test InterventionVector{Int}(2, 2) isa AbstractVector
    @test size(InterventionVector{Int}(2, 2)) == (2, )
    @test size(InterventionVector{Int}(3, 2)) == (3, )
    @test length(InterventionVector{Int}(3, 2)) == 3
end 
@testset "construction without an explicit type" begin
    # test that contructor without a type works, not that it  generates a type `Int`; later
    # tests of `Base.show` demonstrate that the default type is `Int`
    _v1a = InterventionVector(2, 2)
    _v1b = InterventionVector{Int}(2, 2) 
    _v2a = InterventionVector(3, 2)
    _v2b = InterventionVector{Int}(3, 2)
    @test _v1a == _v1b
    @test _v2a == _v2b
end
@testset "offset of zero equal to no offset" begin
    @test InterventionVector(2, 2) == InterventionVector(2, 2; offset=0) 
    @test InterventionVector(3, 2) == InterventionVector(3, 2; offset=0) 
end
@testset "`starttimes` greater that `duration` or equal to 1 equivalent to `nothing`" begin
    @test InterventionVector{Int}(5, 1) == InterventionVector{Int}(5, nothing)
    @test InterventionVector{Int}(5, 6) == InterventionVector{Int}(5, nothing)
    @test InterventionVector{Int}(5, 5; offset=1) == InterventionVector{Int}(5, nothing)
    @test InterventionVector{Int}(5, 2; offset=-1) == InterventionVector{Int}(5, nothing)
    @test InterventionVector{Int}(5, 5) != InterventionVector{Int}(5, nothing)
    @test InterventionVector{Int}(5, 2) != InterventionVector{Int}(5, nothing)
end
@testset "retain offset and raw data in subsequent matrices" begin
    _IV = InterventionVector(10, 1)
    _OIV1 = InterventionVector(_IV; offset=0)
    _OIV2 = InterventionVector(_OIV1; offset=-1)
    _OIV3 = InterventionVector(_OIV2; offset=1)
    _OIV4 = InterventionVector(_OIV3; offset=20)
    _OIV5 = InterventionVector(_OIV4; offset=-15)
    @test _OIV1 == _IV
    @test _OIV2 == InterventionVector(10, nothing)
    @test _OIV3 == _IV
    @test _OIV4 == InterventionVector(10, nothing)
    @test _OIV5 == InterventionVector(10, 6)
end
@testset "construction with explicit non-`Int` type" begin
    @test InterventionVector{Bool}(3, 2)[3]  
    @test !InterventionVector{Bool}(3, 2)[1]  
    @test InterventionVector{Bool}(3, 3; offset=-1)[2]  
end
@testset "invalid arguments to contructor" begin
    @test_throws InexactError InterventionVector(3.3, 2) 
    @test_throws InexactError InterventionVector(3, 2.3)
    @test_throws InexactError InterventionVector(3, 2; offset=0.5)
    @test_throws MethodError InterventionVector(3, [2, 3])
end
@testset "indexing InterventionMatrix" begin
    @test v1[1] == 0
    @test v1[5] == 1
    @test v1[10] == 1
    @test_throws BoundsError v1[0]
    @test_throws BoundsError v1[11]
end
@testset "indexing InterventionMatrix when `starttimes` are nothing" begin
    @test v2[1] == 0
    @test v2[5] == 0
    @test v2[10] == 0
    @test v2[15] == 0
    @test_throws BoundsError v2[0]
    @test_throws BoundsError v2[16]
end
@testset "`_duration`" begin
    @test RenewalDiD._duration(v1) == 10
    @test RenewalDiD._duration(v2) == 15
end
@testset "`_showtimes`" begin
    @test RenewalDiD._showtimes(v1, nothing) == [1, 5, 10]
    @test RenewalDiD._showtimes(v2, nothing) == [1, 15]
end
@testset "output from `Base.show`" begin
    @test repr("text/plain", v1) == v1expectedshowoutput
    @test repr("text/plain", v2) == v2expectedshowoutput
end
