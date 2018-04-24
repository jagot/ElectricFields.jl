function make_field(field_params::Dict{Symbol,Number})
    field_params
end

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
