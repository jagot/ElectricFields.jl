module ElectricFields

using LinearAlgebra
using StaticArrays

using Unitful
using UnitfulAtomic

using AbstractFFTs
using FFTW
using ForwardDiff
using Optim
using Roots
using SpecialFunctions

using IntervalSets

using Parameters

import Base: show
using Format

abstract type AbstractField end

include("units.jl")
include("rotations.jl")
include("relation_dsl.jl")

include("spectra.jl")

include("field_types.jl")
include("time_axis.jl")
include("carriers.jl")
include("envelopes.jl")
include("arithmetic.jl")

include("field_dsl.jl")

# include("sellmeier.jl")
# include("dispersed_fields.jl")

include("strong_field_properties.jl")

end # module
