module ElectricFields

using Unitful
using FFTW
using Parameters
import Base: show
using Printf

include("units.jl")
include("relation_dsl.jl")

include("field_types.jl")
include("time_axis.jl")
include("carriers.jl")
include("envelopes.jl")
include("arithmetic.jl")

include("field_dsl.jl")

include("sellmeier.jl")
include("dispersed_fields.jl")

end # module
