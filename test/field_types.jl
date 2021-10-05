@testset "Field types" begin
    @testset "Linear fields" begin
        @field(F) do
            I₀ = 1.0
            T = 2.0
            σ = 3.0
            Tmax = 3.0
        end

        @test parent(F) === F

        @test field_envelope(F, 0.0) ≈ 1 rtol=1e-7
        @test fluence(F) ≈ 5.452488660809239*3.0*√(2π)/π

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(F) == """
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm"""
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(F) == """
                               Linearly polarized field with
                                 - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
                                   - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
                                   - A₀ = 0.3183 au
                                 – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm"""
        end

        @test ElectricFields.rotation_matrix(F) == Matrix(I, 3, 3)

        t = [0.4, 0.5]
        for fun in (vector_potential, field_amplitude, intensity)
            @test fun(F, t) == [fun(F, t[1]), fun(F, t[2])]
        end
    end

    @testset "Transverse fields" begin
        @field(F) do
            I₀ = 2.0
            T = 2.0
            σ = 3.0
            Tmax = 3.0
            rotation = π/4, [0,0,1]
        end

        @test parent(F) === F

        @test field_envelope(F, 0.0) ≈ √2 rtol=1e-7
        @test fluence(F) ≈ 5.452488660809239*2*3.0*√(2π)/π

        @test intensity(F) == 2.0
        @test amplitude(F) ≈ √2 rtol=1e-7

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(F) == """
                               Transversely polarized field with
                                 - I₀ = 2.0000e+00 au = 7.0188904e16 W cm⁻² =>
                                   - E₀ = 1.4142e+00 au = 727.2178 GV m⁻¹
                                   - A₀ = 0.4502 au
                                 – a LinearTransverseCarrier: Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – and a rotation of 0.25π about [0.000, 0.000, 1.000]
                                 – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm"""
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(F) == """
                               Transversely polarized field with
                                 - I₀ = 2.0000e+00 au = 7.0188904e16 W cm^-2 =>
                                   - E₀ = 1.4142e+00 au = 727.2178 GV m^-1
                                   - A₀ = 0.4502 au
                                 – a LinearTransverseCarrier: Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                 – a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                 – and a rotation of 0.25π about [0.000, 0.000, 1.000]
                                 – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm"""
        end

        @test ElectricFields.rotation_matrix(F) ≈ [1/√2 -1/√2 0
                                                   1/√2 1/√2  0
                                                   0    0     1]

        t = [0.4, 0.5]
        for fun in (vector_potential, field_amplitude, intensity)
            @test fun(F, t) == [fun(F, 0.4)'; fun(F, 0.5)']
        end

        @field(F2) do
            I₀ = 2.0
            T = 2.0
            σ = 3.0
            Tmax = 3.0
            ξ = 1.0
        end

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(F2) == """
                                Transversely polarized field with
                                  - I₀ = 2.0000e+00 au = 7.0188904e16 W cm⁻² =>
                                    - E₀ = 1.4142e+00 au = 727.2178 GV m⁻¹
                                    - A₀ = 0.4502 au
                                  – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                  – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm"""
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(F2) == """
                                Transversely polarized field with
                                  - I₀ = 2.0000e+00 au = 7.0188904e16 W cm^-2 =>
                                    - E₀ = 1.4142e+00 au = 727.2178 GV m^-1
                                    - A₀ = 0.4502 au
                                  – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV)
                                  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
                                  – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm"""
        end
    end

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

    @testset "Ramps" begin
        @field(F) do
            I₀ = 4.0
            tmax = 4.0
            kind = :linear_ramp
        end

       @field(F2) do
            I₀ = 4.0
            tmax = 4.0
            kind = :linear_ramp
            ramp = :down
        end

        withenv("UNITFUL_FANCY_EXPONENTS" => true) do
            @test string(F) == """
                               Linear up-ramp of
                                 - 4.0000 jiffies = 96.7554 as duration, and
                                 - E₀ = 2.0000e+00 au = 1.0284 TV m⁻¹"""
            @test string(F2) == """
                                Linear down-ramp of
                                  - 4.0000 jiffies = 96.7554 as duration, and
                                  - E₀ = 2.0000e+00 au = 1.0284 TV m⁻¹"""
        end
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test string(F) == """
                               Linear up-ramp of
                                 - 4.0000 jiffies = 96.7554 as duration, and
                                 - E₀ = 2.0000e+00 au = 1.0284 TV m^-1"""

            @test string(F2) == """
                                Linear down-ramp of
                                  - 4.0000 jiffies = 96.7554 as duration, and
                                  - E₀ = 2.0000e+00 au = 1.0284 TV m^-1"""
        end

        @test F isa ElectricFields.Ramp

        @test field_amplitude(F, -1) == 0
        @test field_amplitude(F, 0) ≈ 0 rtol=1e-7
        @test field_amplitude(F, 2) ≈ 1 rtol=1e-7
        @test field_amplitude(F, 4-2eps()) ≈ 2 rtol=1e-7
        @test field_amplitude(F, 5) == 0

        @test field_amplitude(F2, -1) == 0
        @test field_amplitude(F2, 0) ≈ 2 rtol=1e-7
        @test field_amplitude(F2, 2) ≈ 1 rtol=1e-7
        @test field_amplitude(F2, 4-2eps()) ≈ 0 atol=1e-10
        @test field_amplitude(F2, 5) == 0

        @test intensity(F, -1) == 0
        @test intensity(F, 0) ≈ 0 rtol=1e-7
        @test intensity(F, 2) ≈ 1 rtol=1e-7
        @test intensity(F, 4-2eps()) ≈ 4 rtol=1e-7
        @test intensity(F, 5) == 0

        @test intensity(F2, -1) == 0
        @test intensity(F2, 0) ≈ 4 rtol=1e-7
        @test intensity(F2, 2) ≈ 1 rtol=1e-7
        @test intensity(F2, 4-2eps()) ≈ 0 atol=1e-10
        @test intensity(F2, 5) == 0

        @test field_amplitude(F, 0, 4.0) ≈ 4.0 rtol=1e-7

        @test polarization(F) == LinearPolarization()

        # @test duration(F) == 4
        @test span(F) == 0..4

        @test isone(period(F))
        @test isone(max_frequency(F))
        @test photon_energy(F) == 2π

        @test dimensions(F) == 1

        @testset "Ramp kind = $(kind)" for (kind,half) in ((:linear_ramp, 0.5),
                                                           (:parabolic_ramp, 0.75),
                                                           (:sin²_ramp, 0.5))
            @field(R) do
                E₀ = 1.0
                tmax = 4.0u"fs"
                kind = kind
            end

            @field(R2) do
                E₀ = 1.0
                tmax = 4.0u"fs"
                kind = kind
                ramp = :down
            end

            @field(C) do
                E₀ = 1.0
                tmax = 4.0u"fs"
                kind = :constant
            end

            F = R + delay(C, duration(R))

            s = span(F)

            @test iszero(field_amplitude(F, s.left-1))
            @test iszero(field_amplitude(F, s.right+1))

            @test field_amplitude(R, austrip(2u"fs")) ≈ half
            @test field_amplitude(R2, austrip(2u"fs")) ≈ half

            @test vector_potential(F, s.right+1) ≈ vector_potential(F, s.right) rtol=1e-7

            t = range(0, stop=1.1span(F).right, length=1000)

            Fv = field_amplitude(F, t)
            Av = vector_potential(F, t)
            Fv2 = -ElectricFields.complex_derivative.(Ref(Base.Fix1(vector_potential, F)), t)

            @test Fv ≈ Fv2 rtol=1e-14
        end
    end
end
