using ElectricFields
using Test
using Unitful
using FFTW

@testset "ElectricFields.jl" begin
    include("namespace_macro.jl")
    include("field_creation.jl")
    include("sellmeier.jl")
end
