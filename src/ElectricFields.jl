module ElectricFields
using Unitful

abstract type AbstractField end

include("parse_block.jl")
include("units.jl")
include("field_types.jl")
include("make_field.jl")

macro field(spec, var)
    spec.head == :-> ||
        error("Expected a block with parameters for definition of the field")
    block = spec.args[2]
    block.head == :block ||
        error("Expected a block with parameters for definition of the field")

    field_params = parse_block(block, Number)
    quote
        $(esc(var)) = make_field($field_params)
    end
end

export @field

end # module
