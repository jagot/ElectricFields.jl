@testset "Constant fields" begin
    @field(F) do
        I₀ = 4.0
        tmax = 4.0
        kind = :constant
    end

    withenv("UNITFUL_FANCY_EXPONENTS" => true) do
        @test string(F) == """
                           Constant field of
                             - 4.0000 jiffies = 96.7554 as duration, and
                             - E₀ = 2.0000e+00 au = 1.0284 TV m⁻¹"""
    end
    withenv("UNITFUL_FANCY_EXPONENTS" => false) do
        @test string(F) == """
                           Constant field of
                             - 4.0000 jiffies = 96.7554 as duration, and
                             - E₀ = 2.0000e+00 au = 1.0284 TV m^-1"""
    end

    @test F isa ElectricFields.ConstantField

    @test field_amplitude(F, -1) == 0
    @test field_amplitude(F, 2) ≈ 2 rtol=1e-7
    @test field_amplitude(F, 5) == 0

    @test intensity(F, -1) == 0
    @test intensity(F, 2) ≈ 4 rtol=1e-7
    @test intensity(F, 5) == 0

    @test field_amplitude(F, 0, 4.0) ≈ 8.0 rtol=1e-7

    @test intensity(F) ≈ 4.0 rtol=1e-7
    @test amplitude(F) ≈ 2.0 rtol=1e-7

    @test polarization(F) == LinearPolarization()

    # @test duration(F) == 4
    @test span(F) == 0..4

    @test isone(period(F))
    @test isone(max_frequency(F))
    @test photon_energy(F) == 2π

    @test dimensions(F) == 1
end
