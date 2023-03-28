using ElectricFields
using Test
using Unitful
using UnitfulAtomic
using IntervalSets
using FFTW
using PrettyTables

function test_approx_eq(a, b; on_fail::Union{Nothing,Function}=nothing, isbroken=false, kwargs...)
    size(a) == size(b) || throw(DimensionMismatch("Cannot compare objects of sizes $(size(a)) and $(size(b))"))
    if !isapprox(a, b; kwargs...)
        @error "Approximate equality failed:"
        na = norm(a)
        nb = norm(b)
        Δ = norm(a-b)
        relΔ = Δ/max(na,nb)
        pretty_table(["|a|" na
                      "|b|" nb
                      "Abs. Δ" Δ
                      "Rel. Δ" relΔ],
                     header=["Quantity", "Value"],
                     hlines=[1], vlines=[],
                     alignment=[:r,:l])
        isnothing(on_fail) || on_fail()
    end

    if !isbroken
        @test isapprox(a, b; kwargs...)
    else
        @test_broken isapprox(a, b; kwargs...)
    end
end

@testset "ElectricFields.jl" begin
    include("namespace_macro.jl")
    include("units.jl")
    include("rotations.jl")
    include("field_creation.jl")
    include("gaussian_pulses.jl")
    include("envelopes.jl")
    include("carriers.jl")
    include("field_types.jl")
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

    include("analytic_continuation.jl")

    include("spectra.jl")

    # include("sellmeier.jl")
end
