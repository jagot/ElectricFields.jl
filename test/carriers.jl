@testset "Carriers" begin
    @field(F1) do
        I₀ = 1.0
        T = 2.0
        τ = 3.0
        σmax = 6.0
        ϕ = π/4
    end

    car1 = carrier(F1)
    @test car1 isa ElectricFields.FixedCarrier
    @test car1(0.4) ≈ sin(π*(0.4+1/4))

    @test string(car1) == "Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV); CEP = 0.25π"

    @test wavelength(car1) |> u"nm" == 2*austrip(1u"c")*u"bohr" |> u"nm"
    @test austrip(period(car1)) == 2.0
    @test austrip(frequency(car1)) == 0.5
    @test wavenumber(car1) |> u"cm^-1" == 1/(2*austrip(1u"c")*u"bohr") |> u"cm^-1"
    @test fundamental(car1) ≈ π
    @test photon_energy(car1) ≈ π
    @test phase(car1) == π/4

    @field(F2) do
        I₀ = 1.0
        T = 2.0
        τ = 3.0
        σmax = 6.0
        ϕ = π/4
        kind = :transverse
    end

    car2 = carrier(F2)
    @test car2 isa ElectricFields.LinearTransverseCarrier

    @test string(car2) == "LinearTransverseCarrier: Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV); CEP = 0.25π"

    @test car2(0.4) ≈ [0,0,sin(π*(0.4+1/4))]
    @test phase_shift(car2, -π/4)(0.4) ≈ [0,0,sin(π*0.4)]

    @test wavelength(car2) |> u"nm" == 2*austrip(1u"c")*u"bohr" |> u"nm"
    @test austrip(period(car2)) == 2.0
    @test austrip(frequency(car2)) == 0.5
    @test wavenumber(car2) |> u"cm^-1" == 1/(2*austrip(1u"c")*u"bohr") |> u"cm^-1"
    @test fundamental(car2) ≈ π
    @test photon_energy(car2) ≈ π
    @test phase(car2) == π/4

    for (ξ,label) in [(0, "0.00"),
                      (1, "1.00 (RCP)"),
                      (-1, "-1.00 (LCP)"),
                      (1/2, "0.50 (right)"),
                      (-1/2, "-0.50 (left)")]
        @field(F3) do
            I₀ = 1.0
            T = 2.0
            τ = 3.0
            σmax = 6.0
            ϕ = π/4
            ξ = ξ
        end
        car3 = carrier(F3)

        @test string(car3) == "Elliptical carrier with ξ = $(label) @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV); CEP = 0.25π"

        @test wavelength(car3) |> u"nm" == 2*austrip(1u"c")*u"bohr" |> u"nm"
        @test austrip(period(car3)) == 2.0
        @test austrip(frequency(car3)) == 0.5
        @test wavenumber(car3) |> u"cm^-1" == 1/(2*austrip(1u"c")*u"bohr") |> u"cm^-1"
        @test fundamental(car3) ≈ π
        @test photon_energy(car3) ≈ π
        @test phase(car3) == π/4
    end
end
