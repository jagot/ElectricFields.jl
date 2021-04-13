@testset "Field arithmetic" begin
    @field(A) do
        I₀ = 1.0
        T = 2.0
        σ = 3.0
        Tmax = 3.0
    end

    @field(A2) do
        I₀ = 1.0
        T = 2.0
        σ = 3.0
        Tmax = 3.0
        ξ = 0.6
    end

    @testset "Summed fields" begin
        @field(B) do
            I₀ = 0.5
            T = 1.0
            σ = 3.0
            Tmax = 3.0
        end

        C = A + B
        @test span(C) == span(A)
        @test span(B) ⊆ span(C)
        @test steps(C, 100) == steps(A, 100/austrip(period(B)))

        @test polarization(C) == LinearPolarization()
        @test dimensions(C) == 1

        @test vector_potential(C, 0.4) == vector_potential(A, 0.4) + vector_potential(B, 0.4)
        @test field_amplitude(C, 0.4) == field_amplitude(A, 0.4) + field_amplitude(B, 0.4)

        # These functions only make sense when adding two field of the same wavelength
        for fun in [wavelength, period, frequency, wavenumber, fundamental, photon_energy]
            @test fun(A+A) == fun(A)
        end

        @test max_frequency(C) == frequency(B)

        @test continuity(C) == Inf

        @test field_amplitude(phase_shift(C, π/3), 0.4) ≈ field_amplitude(phase_shift(A, π/3), 0.4) + field_amplitude(phase_shift(B, π/3), 0.4)

        @test field_amplitude(-A, 0.4) == -field_amplitude(A, 0.4)
        @test field_amplitude(A-B, 0.4) == field_amplitude(A, 0.4) - field_amplitude(B, 0.4)
    end

    @testset "Delayed fields" begin
        B = delay(A, 0.4)
        @test field_amplitude(delay(A,-0.4), 0.0) == field_amplitude(B, 0.8) == field_amplitude(A, 0.4)
        @test field_amplitude(delay(A, π*u"rad"), 0.0) == field_amplitude(A, -1.0)

        @test delay(B) == 0.4
        @test iszero(delay(A))

        @test span(B) == -5.6..6.4
    end

    @testset "Padded fields" begin
        @test_throws ArgumentError PaddedField(A, -3, 1)
        @test_throws ArgumentError PaddedField(A, 3, -1)

        B = PaddedField(A, 3.0u"fs", 8.0u"fs")
        B2 = PaddedField(A2, 3.0u"fs", 8.0u"fs")

        @test parent(B) == A
        @test parent(B2) == A2

        for (F1,F2) in [(A,B), (A2,B2)]
            @test field_amplitude(F2, 0.4) == field_amplitude(F1, 0.4)
            @test field_amplitude(F2, 7.0) == zero(field_amplitude(F1, 0.4))
        end

        @test all(endpoints(span(B)) .≈ (-6-austrip(3.0u"fs"), 6+austrip(8.0u"fs")))

        @test time_integral(B) == time_integral(A)
    end

    @testset "Windowed fields" begin
        B = WindowedField(A, -0.1u"fs", 0.05u"fs")
        B2 = WindowedField(A2, -0.1u"fs", 0.05u"fs")

        @test parent(B) == A
        @test parent(B2) == A2

        for fun in [vector_potential, field_amplitude, intensity]
            for (F1,F2) in [(A,B), (A2,B2)]
                @test field_amplitude(F2, 0.4) == field_amplitude(F1, 0.4)
                @test field_amplitude(F2, -5.0) == zero(field_amplitude(F1, 0.4))
                @test field_amplitude(F2, 3.0) == zero(field_amplitude(F1, 0.4))
            end
        end

        @test all(endpoints(span(B)) .≈ (-austrip(0.1u"fs"), austrip(0.05u"fs")))

        @test field_amplitude(phase_shift(B, 2.0), 0.4) == field_amplitude(phase_shift(A, 2.0), 0.4)
        @test field_amplitude(phase_shift(B, 2.0), 3.0) == zero(field_amplitude(phase_shift(A, 2.0), 0.4))
    end
end
