module ElectricFields
using Unitful

abstract type AbstractField end

include("parse_block.jl")
include("units.jl")
include("field.jl")

end # module
