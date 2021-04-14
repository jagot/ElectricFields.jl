import ElectricFields: Iaustrip

@testset "Field creation" begin
    @testset "Parametric field creation" begin
        #    Test parametric creation of fields, both from inside a function and
        #    from global scope.

        function test_parametric_fields(l)
            Is = range(0u"W/cm^2",stop=1e6u"W/cm^2", length=l)

            for (i,I₀) in enumerate(Is)
                @field(IR) do
                    λ = 800.0u"nm"
                    I₀ = I₀
                    τ = 6.2u"fs"
                    Tmax = 5
                end
                @test intensity(IR) == Iaustrip(I₀)
            end
        end
        test_parametric_fields(11)
        test_parametric_fields(8)

        Is = range(0u"W/cm^2",stop=1e6u"W/cm^2", length=11)

        for (i,I₀) in enumerate(Is)
            @field(IR) do
                λ = 800.0u"nm"
                I₀ = I₀
                τ = 6.2u"fs"
                Tmax = 5
            end
            @test intensity(IR) == Iaustrip(I₀)
        end

        @field(IR) do
            λ = 800.0u"nm"
            I₀ = 1e14u"W/cm^2"
            τ = 6.2u"fs"
            Tmax = 5
        end

        @test dimensions(IR) == 1

        @field(IR2) do
            λ = 800.0u"nm"
            I₀ = 1e14u"W/cm^2"
            τ = 6.2u"fs"
            Tmax = 5
            ξ = 0.4
        end

        @test dimensions(IR2) == 3
    end

    @test_throws ArgumentError @field(F) do
        I₀ = 1e14u"W/cm^2"
        λ = 800.0u"nm"
        kind = :there_is_no_such_field
    end

    @testset "Competing quantities" begin
        @field(A) do
            I₀ = 1e14u"W/cm^2"
            λ = 800.0u"nm"
            τ = 6.2u"fs"
            σmax = 4
        end

        @test period(A) |> u"fs" ≈ 2.6685127615852164u"fs"
        @test wavenumber(A) |> u"cm^-1" ≈ 12500.0u"cm^-1"
        @test frequency(A) |> u"THz" ≈ 374.7405725u"THz"
        @test photon_energy(A) ≈ 0.05695419066118882

        @test ponderomotive_potential(A) |> u"eV" ≈ 5.9758653162292505u"eV"
        @test intensity(A) ≈ 0.0028494532412131697
        @test amplitude(A) ≈ 0.05338026769639452
        @test vector_potential(A) ≈ 0.9372491659822719

        @field(B) do
            I₀ = 1e14u"W/cm^2"
            T = 2.6685127615852164u"fs"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(C) do
            I₀ = 1e14u"W/cm^2"
            ν = 12500.0u"cm^-1"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(D) do
            I₀ = 1e14u"W/cm^2"
            f = 374.7405725u"THz"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(E) do
            I₀ = 1e14u"W/cm^2"
            ω = 2.354564459136067u"Prad/s"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(E) do
            I₀ = 1e14u"W/cm^2"
            ħω = 0.05695419066118882
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(F) do
            I₀ = 0.0028494532412131697
            λ = 800.0u"nm"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(G) do
            E₀ = 0.05338026769639452
            λ = 800.0u"nm"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(H) do
            A₀ = 0.9372491659822719
            λ = 800.0u"nm"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(I) do
            Uₚ = 5.9758653162292505u"eV"
            λ = 800.0u"nm"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(J) do
            Uₚ = 0.21960899978361606
            λ = 800.0u"nm"
            τ = 6.2u"fs"
            σmax = 4
        end

        @field(K) do
            I₀ = 1e14u"W/cm^2"
            λ = 800.0u"nm"
            τ = 6.2u"fs"
            σmax = 4
            ξ = 1.0
        end

        @testset "$(fun)" for (fun,u) in [(wavelength, u"nm"),
                                          (period, u"fs"),
                                          (frequency, u"THz"),
                                          (photon_energy, NoUnits),
                                          (ponderomotive_potential, u"eV"),
                                          (intensity, NoUnits),
                                          (amplitude, NoUnits),
                                          (vector_potential, NoUnits)]
            @testset "i = $i" for (i,f) in enumerate([B,C,D,E,F,G,H,I,J,K])
                @test fun(f) |> u ≈ fun(A) |> u
            end
        end
    end
end
