using ElectricFields
using Test
using Unitful
using UnitfulAtomic
using IntervalSets
using FFTW

@testset "ElectricFields.jl" begin
    include("namespace_macro.jl")
    include("units.jl")
    include("rotations.jl")
    include("field_creation.jl")
    include("gaussian_pulses.jl")
    include("envelopes.jl")
    include("carriers.jl")
    include("constant_fields.jl")
    include("arithmetic.jl")

    @testset "Time axis" begin
        @field(F) do
            I₀ = 1.0
            T = 2.0
            σ = 3.0
            Tmax = 3.0
        end

        @test ElectricFields.default_sampling_frequency(F) == 100/2

        @test steps(F) == 600
        @test steps(F, 300) == 1800

        @test timeaxis(F) == range(-6, stop=6, length=600)
        @test timeaxis(F, 300) == range(-6, stop=6, length=1800)
        @test timeaxis(F, 0) == -6:12:-6
        @test timeaxis(F, (26,0.5)) == range(-6, length=26, step=0.5)
        @test timeaxis(F, (24,0.5)) == range(-6, length=24, step=0.5)
    end

    @testset "Strong-field properties" begin
        @field(F) do
            I₀ = 1.0
            T = 2.0
            σ = 3.0
            Tmax = 3.0
        end

        @test austrip(ponderomotive_potential(F)) ≈ 1/(4*π^2) rtol=1e-7
        @test austrip(ponderomotive_potential(F+F)) ≈ 2/(4*π^2) rtol=1e-7

        @test keldysh(F, 2u"Eh_au") ≈ 2π rtol=1e-7

        @test free_oscillation_amplitude(F) ≈ 1/π^2 rtol=1e-7
        @test free_oscillation_amplitude(delay(F, 100.0)) ≈ 1/π^2 rtol=1e-7
        @test free_oscillation_amplitude(F+F) ≈ 2/π^2 rtol=1e-7
    end

    # include("sellmeier.jl")
end
