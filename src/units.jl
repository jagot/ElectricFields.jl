@derived_dimension ElectricField Unitful.𝐈^-1*Unitful.𝐋*Unitful.𝐌*Unitful.𝐓^-3
@derived_dimension Intensity Unitful.𝐌*Unitful.𝐓^-3

# Unit arithmetic is performed "outside" of @u_str, to ensure type
# stability and allow precompilation
base_units = Dict{Symbol,Unitful.FreeUnits}(
    :λ => u"nm",
    :I₀ => u"W"/u"cm"^2,
    :E₀ => u"V"/u"m",
    :τ => u"fs",
    :σ => u"fs",
    :σ′ => u"fs",
    :σmax => NoUnits,
    :σ′max => NoUnits,
    :toff => u"fs",
    :tmax => u"fs",
    :Tmax => NoUnits,
    :T => u"fs",
    :f => u"THz",
    :ν => u"cm"^-1,
    :ω => u"Trad"/u"s",
    :ħω => u"eV",
    :Uₚ => u"eV"
)

function get_unitful_quantity(field_params::Dict{Symbol,Any}, sym::Symbol)
    v = field_params[sym]
    typeof(v) <: Quantity || sym ∉ keys(base_units) ? v : v*base_units[sym]
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

usplit(u) = ustrip(u),unit(u)
