# test functions related to generating and using `InterventionMatrix` structs

using RenewalDiD
using Test

tw1 = "InterventionArray with no intervention in any group"
tw2 = "All groups in InterventionArray have intervention before end of duration"

@testset "constructing an InterventionArray" begin
    @test InterventionArray{Int}(2, [[2, nothing], [1, 2]]) isa AbstractInterventionArray
    # (behaviour of `mutewarnings` is tested below)
    @test InterventionArray{Int}(2, [[2, nothing], [1, 2]]) isa AbstractArray
    @test size(InterventionArray{Int}(2, [[2, nothing], [1, 2]])) == (2, 2, 2)
    @test size(InterventionArray{Int}(3, [[2, nothing], [1, 2]])) == (3, 2, 2)
    @test size(InterventionArray{Int}(3, [[2, nothing, 2], [1, 2, 2]])) == (3, 3, 2)
    @test size(InterventionArray{Int}(3, [[2, nothing], [1, 2], [nothing, 3]])) == (3, 2, 3)
end
@testset "construction without an explicit type" begin
    # test that contructor without a type works, not that it  generates a type `Int`; later
    # tests of `Base.show` demonstrate that the default type is `Int`
    _M1a = InterventionArray(2, [[2, nothing], [1, 2]])
    _M1b = InterventionArray(3, [[2, nothing, 2], [1, 2, 2]])
    _M2a = InterventionArray{Int}(2, [[2, nothing], [1, 2]])
    _M2b = InterventionArray{Int}(3, [[2, nothing, 2], [1, 2, 2]])
    @test _M1a == _M2a
    @test _M1b == _M2b
end
@testset "construction with different arrangement of arguments" begin
    _M1a = InterventionArray(2, [[2, nothing], [1, 2]])
    _M1b = InterventionArray(3, [[2, nothing, 2], [1, 2, 2]])
    _M2a = InterventionArray(2, [2, nothing], [1, 2])
    _M2b = InterventionArray(3, [2, nothing, 2], [1, 2, 2])
    _M3a = InterventionArray(2, [2  1; nothing  2])
    _M3b = InterventionArray(3, [2  1; nothing  2; 2  2])
    @test _M1a == _M2a
    @test _M1a == _M3a
    @test _M1b == _M2b
    @test _M1b == _M3b
end
@testset "retain offset and raw data in subsequent matrices" begin
    _IM = InterventionArray(10, [1, 5, 9], [11, 15, nothing])
    _OM1 = InterventionArray(_IM; offset=0)
    _OM2 = InterventionArray(_OM1; offset=0)
    _OM3 = InterventionArray(_OM2; offset=-1)
    _OM4 = InterventionArray(_OM3; offset=1)
    _OM5 = InterventionArray(_OM4; offset=1)
    _OM6 = InterventionArray(_OM5; offset=-9)
    @test _IM == _OM1
    @test _OM1 == _OM2
    @test _OM3 == InterventionArray(10, [nothing, 4, 8], [10, nothing, nothing])
    @test _OM4 == _IM
    @test _OM5 == InterventionArray(10, [2, 6, 10], [nothing, nothing, nothing])
    @test _OM6 == InterventionArray(10, [nothing, nothing, 1], [3, 7, nothing])
end 
@testset "create by only changing offsets" begin
    A1 = InterventionArray(10, [1, 5, 9]; offset=[0, -10, -5, 5, 10])
    A2 = InterventionArray(
        10,
        [nothing, 5, 9],
        [nothing, nothing, nothing],
        [nothing, nothing, 4],
        [6, 10, nothing],
        [nothing, nothing, nothing]
    )
    @test A1 == A2
end
@testset "error if interventions and offsets do not match" begin
    @test_throws ArgumentError InterventionArray(
        10, [1, 5, 9], [3, 5, 8]; 
        offset=[0, -10, 10]
    )
end
@testset "warnings during construction" begin
    @test_warn tw1 InterventionArray(2, [nothing, nothing], [3, 3])
    @test_nowarn InterventionArray(2, [nothing, nothing], [3, 3]; mutewarnings=true)
    @test_warn tw1 InterventionArray(2, [nothing, nothing], [3, 3]; mutewarnings=false)
    @test_warn tw2 InterventionArray(5, [2, 3], [4, 5])
    @test_nowarn InterventionArray(5, [2, 3], [4, 5]; mutewarnings=true)
    @test_warn tw2 InterventionArray(5, [2, 3], [4, 5]; mutewarnings=false)
    @test_nowarn InterventionArray(2, [2, 3], [2, 2])
end 
@testset "construction with explicit non-`Int` type" begin
    @test InterventionArray{Bool}(3, [2, 4], [1, 1])[3, 1, 1]  
    @test !InterventionArray{Bool}(3, [2, 4], [1, 1])[3, 1, 2]  
    @test InterventionArray{Bool}(3, [2, 4], [1, 1]; offset=1)[3, 1, 2]  
end
@testset "invalid arguments to contructor" begin
    @test_throws InexactError InterventionArray(3.3, [2, 4], [2, 4]) 
    @test_throws InexactError InterventionArray(3, [2.3, 2], [nothing, nothing])
    @test_throws InexactError InterventionArray(3, [2, 4], [2, 4]; offset=0.5)
    @test_throws InexactError InterventionArray(3, [2, 4], [2, 4]; offset=[0, 0.5])
    @test_throws MethodError InterventionArray(3, [[2, 4], [2  3; 3  5]])
    @test_throws MethodError InterventionArray(3, [2, 4, [2, 3], [3, 5]])
    @test_throws MethodError InterventionArray(3, 2, 4, [2, 3], [3, 5])
end
@testset "indexing InterventionMatrix" begin
    A1 = InterventionArray{Int}(5, [2, 3], [1, 3], [3, nothing]) 
    @test A1[1] == 0
    @test A1[1, 1, 1] == 0
    @test A1[3, 1, 1] == 1
    @test A1[1, 2, 1] == 0
    @test_throws BoundsError A1[0]
    @test_throws BoundsError A1[0, 1, 1]
    @test_throws BoundsError A1[1, 0, 1]
    @test_throws BoundsError A1[1, 1, 0]
    @test_throws BoundsError A1[6, 1, 1]
    @test_throws BoundsError A1[1, 3, 1]
    @test_throws BoundsError A1[1, 1, 4]
end
@testset "`_duration`" begin
    @test RenewalDiD._duration(InterventionArray{Int}(5, [2, 3], [1, 3], [3, nothing])) == 5
end
@testset "concatenate matrices and arrays" begin
    v1 = InterventionVector(10, 2)
    A1 = InterventionArray(10, [2, 7], [3, nothing])
    A2 = InterventionArray(10, [2, 7], [3, nothing]; offset=[3, -2])
    A3 = InterventionArray(9, [2, 7], [3, nothing]; offset=[3, -2])
    A4 = InterventionArray(10, [2, 7, 5], [3, nothing, nothing])
    M1 = InterventionMatrix(10, [2, nothing])
    M2 = InterventionMatrix(10, [2, nothing]; offset=1)
    M3 = InterventionMatrix(8, [2, nothing]; offset=1)
    M4 = InterventionMatrix(10, [2, nothing, 2])
    N1 = cat(A1; dims=2)
    N2 = cat(M1; dims=3)
    N2a = cat(M1; dims=Val{3}())
    N2b = interventioncat(M1)
    N3 = cat(A1, M1; dims=3)
    N3a = interventioncat(A1, M1)
    N4 = cat(M1, A1; dims=3)
    N4a = interventioncat(M1, A1)
    N5 = cat(A1, N4; dims=3)
    N5a = interventioncat(A1, N4)
    N6 = cat(M1, M2; dims=3)
    N6a = interventioncat(M1, M2)
    N7 = cat(A1, A2; dims=3)
    N7a = interventioncat(A1, A2)
    @test_throws ArgumentError cat(A1; dims=1)
    @test_throws ArgumentError cat(A1; dims=Val{1}())
    @test_throws ArgumentError cat(A1, A2; dims=Val{1}())
    @test_throws ArgumentError cat(A1, M1; dims=Val{1}())
    @test_throws ArgumentError cat(M1, A1; dims=Val{1}())
    @test_throws MethodError cat(v1, A1; dims=Val{3}())
    @test_throws MethodError cat(v1, M1; dims=Val{3}())
    @test_throws ArgumentError cat(A1, v1; dims=Val{3}())
    @test_throws ArgumentError cat(M1, v1; dims=Val{3}())
    @test_throws ArgumentError interventioncat(M1, v1)
    @test N1 == A1
    @test N2 == InterventionArray(10, [2, nothing])
    @test N2a == InterventionArray(10, [2, nothing])
    @test N2b == InterventionArray(10, [2, nothing])
    @test N3 == InterventionArray(10, [2, 7], [3, nothing], [2, nothing])
    @test N3a == InterventionArray(10, [2, 7], [3, nothing], [2, nothing])
    @test N4 == InterventionArray(10, [2, nothing], [2, 7], [3, nothing])
    @test N4a == InterventionArray(10, [2, nothing], [2, 7], [3, nothing])
    @test N5 == InterventionArray(10, [2, 7], [3, nothing], [2, nothing], [2, 7], [3, nothing])
    @test N5a == InterventionArray(10, [2, 7], [3, nothing], [2, nothing], [2, 7], [3, nothing])
    @test N6 == InterventionArray(10, [2, nothing], [2, nothing]; offset=[0, 1])
    @test N6a == InterventionArray(10, [2, nothing], [2, nothing]; offset=[0, 1])
    @test N6 == InterventionArray(10, [2, nothing], [3, nothing])
    @test N7 == InterventionArray(
        10, [2, 7], [3, nothing], [2, 7], [3, nothing]; 
        offset=[0, 0, 3, -2]
    )
    @test N7a == N7
    @test_throws DimensionMismatch cat(A1, A3; dims=3)
    @test_throws DimensionMismatch interventioncat(A1, A3)
    @test_throws DimensionMismatch cat(A1, A4; dims=3)
    @test_throws DimensionMismatch cat(A1, M3; dims=3)
    @test_throws DimensionMismatch cat(A1, M4; dims=3)
    @test_throws DimensionMismatch cat(A3, A1; dims=3)
    @test_throws DimensionMismatch cat(A4, A1; dims=3)
    @test_throws DimensionMismatch cat(M3, A1; dims=3)
    @test_throws DimensionMismatch cat(M4, A1; dims=3)
end
