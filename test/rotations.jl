using StaticArrays
using LinearAlgebra
import ElectricFields: compute_rotation, rotation_angle, rotation_axis

@testset "Rotations" begin
    @test compute_rotation(nothing) == I
    @test compute_rotation(π*I) == I

    @testset "From axis–angle" begin
        @test compute_rotation((0,[0,0,1])) == I
        @test compute_rotation((π/4,[0,0,1])) ≈ [1/√2 -1/√2 0
                                                 1/√2 1/√2  0
                                                 0    0     1]
        @test compute_rotation((π/2,[0,0,1])) ≈ [0 -1 0
                                                 1  0 0
                                                 0  0 1]
        @test compute_rotation((π,[0,0,1])) ≈ [-1 0 0
                                               0 -1 0
                                               0  0 1]
        @test compute_rotation((π/2,[0,1,1])) ≈ [0    -1/√2 1/√2
                                                 1/√2  1/2  1/2
                                                 -1/√2 1/2  1/2]
    end

    @testset "From rotation matrix" begin
        @test compute_rotation([1.0 1.0 0.0
                                0.0 1.0 1.0
                                0.0 0.0 1.0]) ≈ I rtol=1e-14
        @test compute_rotation([1.0 1.0 1.0
                                0.0 1.0 1.0
                                0.0 0.0 1.0]) ≈ I rtol=1e-14

        @test_throws DimensionMismatch compute_rotation([1.0 0.0
                                                         0.0 1.0])

        @test_throws ArgumentError compute_rotation([1.0 1.0 0.0
                                                     0.0 0.0 0.0
                                                     0.0 0.0 1.0])
    end

    @testset "Rotation axis and angle from matrix" begin
        @test rotation_angle(Matrix(I, 3, 3)) == 0
        for (angle,axis) in [(π/4,[0,0,1]),
                             (π/3,[0,0,1]),
                             (π/2,[0,0,1]),
                             (π,[0,0,1]),
                             (π/2,[0,1,1])]
            R = compute_rotation((angle, axis))
            @test rotation_angle(R) ≈ angle
            @test rotation_axis(R) ≈ normalize(axis)
        end
    end
end
