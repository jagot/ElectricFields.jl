import ElectricFields: time_integral

@testset "Gaussian pulses" begin
    @testset "Equivalent parameters" begin
        @field(F1) do
            I₀ = 1.0
            T = 2.0
            τ = 3.0
            σmax = 6.0
        end
        env1 = envelope(F1)
        @test env1 isa ElectricFields.GaussianEnvelope
        σ1 = 3.0/(2*√(2log(2)))
        @test 6σ1 ≤ env1.tmax ≤ 6σ1 + 2
        @test time_integral(F1) ≈ σ1*√(2π)
        @test string(env1) == "Gaussian envelope of duration 72.5665 as (intensity FWHM; ±6.28σ)"

        @field(F2) do
            I₀ = 1.0
            T = 2.0
            σ = 3.0
            tmax = 6.0
        end
        env2 = envelope(F2)
        @test env2 isa ElectricFields.GaussianEnvelope
        τ2 = 3.0*(2*√(2log(2)))
        @test env2.τ ≈ τ2
        @test 6.0 ≤ env2.tmax ≤ 6.0 + 2.0

        @field(F3) do
            I₀ = 1.0
            T = 2.0
            σ = 3.0
            Tmax = 3.0
        end
        env3 = envelope(F3)
        @test env3 isa ElectricFields.GaussianEnvelope
        τ2 = 3.0*(2*√(2log(2)))
        @test env3.τ ≈ τ2
        @test env3.tmax == 6.0

        @test continuity(F3) == Inf
        @test span(F3) == -6.0..6.0
        @test time_integral(F3) ≈ 3.0*√(2π)
    end

    @testset "FWHM derivation" begin
        # FWHM is defined for the intensity envelope, which is given
        # by the field amplitude envelope squared, which in turn
        # depends on the derivative of the vector amplitude. Here we
        # check that the creation routine for Gaussian pulses manages
        # to find an FWHM for the vector amplitude that is consistent
        # with the requested FWHM for the intensity envelope.
        @field(F) do
            I₀ = 1.0
            T = 2.0
            τ = 3.0
            σmax = 6.0
        end
        @test intensity(F, 1.5) ≈ 0.5 rtol=1e-7

        @field(F2) do
            I₀ = 1.0
            T = 2.0
            τ = 3.0
            toff = 4.0
            tmax = 5.0
            env = :trunc_gauss
        end
        @test intensity(F2, 1.5) ≈ 0.5 rtol=1e-7
    end

    @testset "Smooth turn-off of truncated Gaussian" begin
        @field(F) do
            I₀ = 1.0
            T = 2.0
            τ = 3.0
            toff = 4.0
            tmax = 5.0
            env = :trunc_gauss
        end
        env = envelope(F)
        @test string(env) == "Truncated Gaussian envelope of duration 3.0000 jiffies = 72.5665 as (intensity FWHM; turn-off from 96.7554 as to 120.9442 as)"

        α = env.α
        for t = -3.5:0.5:3.5
            @test env(t) ≈ exp(-α*t^2) rtol=1e-12
        end

        @test span(F) == -5..5
    end
end
