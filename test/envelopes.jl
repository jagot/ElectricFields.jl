@testset "Envelopes" begin
    @testset "Trapezoidal envelopes" begin
        @field(F) do
            I₀ = 1.0
            T = 2.0
            ramp_up = 1.0
            flat = 3.0
            ramp_down = 2.0
            env = :trapezoidal
        end

        env = envelope(F)
        @test env isa ElectricFields.TrapezoidalEnvelope

        @test env(-1.0) == 0.0
        @test env(0.0) == 0.0
        @test env(1.0) == 0.5
        @test env(2.5) == 1.0
        @test env(7.9) == 1.0
        @test env(10.0) == 0.5
        @test env(12.0) == 0.0
        @test env(13.0) == 0.0

        @test continuity(env) == 0
        @test span(env) == 0..12

        @test string(env) == "/1‾3‾2\\ cycles trapezoidal envelope"

        @field(F2) do
            I₀ = 1.0
            T = 2.0
            ramp = 1.0
            flat = 3.0
            env = :trapezoidal
        end

        env2 = envelope(F2)

        @test env2.ramp_up == env2.ramp_down == 1.0
        @test span(env2) == 0..10

        @test_throws ErrorException ElectricFields.TrapezoidalEnvelope(-1.0, 1.0, 1.0, 1.0)
        @test_throws ErrorException ElectricFields.TrapezoidalEnvelope(1.0, -1.0, 1.0, 1.0)
        @test_throws ErrorException ElectricFields.TrapezoidalEnvelope(1.0, 1.0, -1.0, 1.0)
        @test_throws ErrorException ElectricFields.TrapezoidalEnvelope(0.0, 0.0, 0.0, 1.0)
    end

    @testset "cos² envelopes" begin
        @field(F) do
            I₀ = 1.0
            T = 2.0
            cycles = 3.0
            env = :cos²
        end

        env = envelope(F)

        @test env(0.0) == 1.0
        @test env(6/4) ≈ 1/2
        @test env(3.0) == 0.0
        @test env(3.1) == 0.0

        @test string(env) == "3.00 cycles cos² envelope"

        @test continuity(env) == 0
        @test span(env) == -3..3
    end

    @testset "Gaussian envelopes" begin
        for τ ∈ range(0.1, stop=2.0, length=11)
            @field(F) do
                I₀ = 1.0
                T = 1.0
                τ = τ
                σmax = 6.0
            end

            @test intensity(F, τ/2) ≈ 1/2 rtol=1e-5
        end
    end
end
