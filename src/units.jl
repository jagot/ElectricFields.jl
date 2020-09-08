@derived_dimension ElectricField Unitful.𝐈^-1*Unitful.𝐋*Unitful.𝐌*Unitful.𝐓^-3
@derived_dimension VectorPotential Unitful.𝐈^-1*Unitful.𝐋*Unitful.𝐌*Unitful.𝐓^-2
@derived_dimension Intensity Unitful.𝐌*Unitful.𝐓^-3

# Unit arithmetic is performed "outside" of @u_str, to ensure type
# stability and allow precompilation
const base_units = Dict{Symbol,Unitful.FreeUnits}(
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
    typeof(v) <: Quantity || sym ∉ keys(base_units) ? v : v*aunit(base_units[sym])
end

usplit(u) = ustrip(u),unit(u)
