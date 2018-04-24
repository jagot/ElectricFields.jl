using Unitful

function set_base_units(unit_specs::Dict{Symbol, Any})
    println(unit_specs)
end

macro set_base_units(spec)
    spec.head == :-> ||
        error("Expected a block with parameters for definition of the field")
    block = spec.args[2]
    block.head == :block ||
        error("Expected a block with parameters for definition of the field")

    set_base_units(parse_block(block, Any))
end

export @set_base_units
