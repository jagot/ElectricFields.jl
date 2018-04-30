using Unitful

base_units = Dict(:λ => u"nm",
                  :I₀ => u"W/cm^2",
                  :E₀ => u"V/m",
                  :τ => u"fs",
                  :T => u"s",
                  :f => u"Hz",
                  :ν => u"cm^-1",
                  :ω => u"rad/s")

function get_unitful_quantity(field_params::Dict{Symbol,Union{Number,Quantity}}, sym::Symbol)
    v = field_params[sym]
    typeof(v) <: Quantity ? v : v*base_units[sym]
end

function set_base_units(unit_specs::Dict{Symbol, Any})
    global base_units
    unknown_units = setdiff(keys(unit_specs),
                            keys(base_units))
    length(unknown_units) != 0 && error("Unknown base unit, $(join(unknown_units, ", "))")

    for (k,v) in unit_specs
        base_units[k] = v
    end

    nothing
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
