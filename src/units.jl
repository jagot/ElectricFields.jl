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

const Iau = 3.5094452e16u"W"/u"cm"^2

# We have to use this hack, since UnitfulAtomic does not consider Iau
# as the atomic unit of intensity (which in itself is a contested name
# for the quantity in question).
oneaunit(x) = aunit(x)
oneaunit(x::Intensity) = Iau

function get_unitful_quantity(field_params::Dict{Symbol,Any}, sym::Symbol)
    v = field_params[sym]
    typeof(v) <: Quantity || sym ∉ keys(base_units) ? v : v*oneaunit(1base_units[sym])
end

function shift_unit(u::U, d, n) where {U<:Unitful.Unit}
    i₀ = floor(Int, d/n)
    d = 3i₀
    iszero(d) && return u,0
    for (i,tens) in enumerate(Unitful.tens(u) .+ (d:(-3*sign(d)):0))
        haskey(Unitful.prefixdict, tens) && return U(tens, u.power),i-1 + i₀
    end
    u,0
end

function shift_unit(u::Unitful.FreeUnits, d)
    tu = typeof(u)
    us,ds,a = tu.parameters
    ui = findfirst(u -> u.power>0, us)

    uu,i = shift_unit(us[ui], d, 3)

    Unitful.FreeUnits{(us[1:ui-1]...,uu,us[ui+1:end]...), ds, a}(),i
end

function si_round(q::Quantity; fspec="{1:.4f} {2:s}")
    v,u = ustrip(q), unit(q)
    if !iszero(v)
        u,i = shift_unit(u, log10(abs(v)))
        q = u(q)
    end
    format(fspec, ustrip(q), unit(q))
end

au2si(v, u) = u(v*aunit(u))
au2si(u::Unitful.Units)= au2si(1, u)

au2si_round(v, u; kwargs...) =
    si_round(au2si(v, u); kwargs...)

I2si_round(I; kwargs...) = si_round(I*Iau)

Iaustrip(I::Intensity) = NoUnits(I/Iau)
