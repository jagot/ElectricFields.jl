using ElectricFields
using Test
using Unitful
using UnitfulAtomic
using FFTW

@testset "ElectricFields.jl" begin
    include("namespace_macro.jl")
    include("units.jl")
    include("rotations.jl")
    include("field_creation.jl")
    include("gaussian_pulses.jl")
    # include("sellmeier.jl")
end
