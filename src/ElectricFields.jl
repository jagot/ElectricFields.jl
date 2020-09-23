module ElectricFields

using LinearAlgebra
using StaticArrays

using Unitful
using UnitfulAtomic

# using FFTW
using ForwardDiff
using Optim

using IntervalSets

using Parameters

import Base: show
using Formatting

include("units.jl")
include("relation_dsl.jl")

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
