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

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(C) == """
                               ┌ Linearly polarized field with
                               │   - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
                               │     - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
                               │     - A₀ = 0.3183 au
                               │   – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                               │   – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                               │   – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
                               ⊕
                               │ Linearly polarized field with
                               │   - I₀ = 5.0000e-01 au = 1.7547226e16 W cm⁻² =>
                               │     - E₀ = 7.0711e-01 au = 363.6089 GV m⁻¹
                               │     - A₀ = 0.1125 au
                               │   – a Fixed carrier @ λ = 7.2516 nm (T = 24.1888 as, ω = 6.2832 Ha = 170.9742 eV)
                               │   – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±1.00σ)
                               └   – Uₚ = 0.0032 Ha = 86.1591 meV => α = 0.0179 Bohr = 947.8211 fm
                               """
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(C) == """
                               ┌ Linearly polarized field with
                               │   - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
                               │     - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
                               │     - A₀ = 0.3183 au
                               │   – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                               │   – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                               │   – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
                               ⊕
                               │ Linearly polarized field with
                               │   - I₀ = 5.0000e-01 au = 1.7547226e16 W cm^-2 =>
                               │     - E₀ = 7.0711e-01 au = 363.6089 GV m^-1
                               │     - A₀ = 0.1125 au
                               │   – a Fixed carrier @ λ = 7.2516 nm (T = 24.1888 as, ω = 6.2832 Ha = 170.9742 eV)
                               │   – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±1.00σ)
                               └   – Uₚ = 0.0032 Ha = 86.1591 meV => α = 0.0179 Bohr = 947.8211 fm
                               """
        end
    end

    @testset "Delayed fields" begin
        B = delay(A, 0.4)
        @test field_amplitude(delay(A,-0.4), 0.0) == field_amplitude(B, 0.8) == field_amplitude(A, 0.4)
        @test field_amplitude(delay(A, π*u"rad"), 0.0) == field_amplitude(A, -1.0)

        @test delay(B) == 0.4
        @test iszero(delay(A))

        @test span(B) == -5.6..6.4

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(B) == """
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
                                 – delayed by 0.4000 jiffies = 9.6755 as"""
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(B) == """
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
                                 – delayed by 0.4000 jiffies = 9.6755 as"""
        end
    end

    @testset "Combined" begin
        @field(X) do
            I₀ = 1.0
            ω = 1.0
            cycles = 6.0
            env = :cos²
            ξ = 1.0
        end

        @field(Y) do
            I₀ = 1.0
            ω = 2.0
            cycles = 6.0
            env = :cos²
            ξ = -1.0
        end

        Z = X + delay(Y, 3/2π)

        @test field_amplitude(Z, 4.0) ≈ field_amplitude(X, 4.0) + field_amplitude(Y, 4.0 - 3/2π)
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

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(B) == """
                               Padding before 124.0241 jiffies = 3.0000 fs and after 330.7310 jiffies = 8.0000 fs of
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm"""
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(B) == """
                               Padding before 124.0241 jiffies = 3.0000 fs and after 330.7310 jiffies = 8.0000 fs of
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm"""
        end
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

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(B) == """
                               Window from -4.1341 jiffies = -100.0000 as to 2.0671 jiffies = 50.0000 as of
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm"""
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(B) == """
                               Window from -4.1341 jiffies = -100.0000 as to 2.0671 jiffies = 50.0000 as of
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm"""
        end
    end
end
