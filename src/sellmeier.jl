# * Media

using Unitful
import Unitful: Length, Area, Frequency
import Base: show
using Printf

@doc raw"""
    Medium(B, C)

The [Sellmeier
equations](https://en.wikipedia.org/wiki/Sellmeier_equation) are used
to describe dispersion in glass using a series of resonances:

```math
n^2(\lambda) =
1 + \sum_i \frac {B_i\lambda^2}{\lambda^2-C_i}.
```
"""
struct Medium{T<:Real,A<:Area{T}}
    B::Vector{T}
    C::Vector{A}
end

function show(io::IO, m::Medium)
    write(io, "Medium(")
    show(io, m.B)
    write(io, @sprintf(", [%s])",
                       join(string.(m.C), ", ")))
end

# ** BK7
"""
    BK7

[Borosilicate glass](https://en.wikipedia.org/wiki/Borosilicate_glass)
is commonly used in optical lenses
"""
BK7 = Medium([1.03961212, 0.231792344, 1.01046945],
             [6.00069867e-3u"μm"^2, 2.00179144e-2u"μm"^2, 1.03560653e2u"μm"^2])

# ** SiO₂

"""
    SiO₂

[Fused silica](https://en.wikipedia.org/wiki/Fused_quartz#Optical_properties)
"""
SiO₂ = Medium([0.6961663, 0.4079426, 0.8974794],
              [0.0684043u"μm", 0.1162414u"μm", 9.896161u"μm"].^2)

# * Refractive index

function n²(m::Medium, λ::Length)
    if !isfinite(λ)
        1.0
    else
        1.0 + sum([Bᵢ*λ^2/(λ^2-Cᵢ)
                   for (Bᵢ,Cᵢ) ∈ zip(m.B, m.C)]) |> NoUnits
    end
end

maybe_complex(x) = x<0 ? complex(x) : x
refractive_index(m::Medium, λ::Length) = √(maybe_complex(n²(m, λ)))
refractive_index(m::Medium, f::Frequency) = refractive_index(m, u"c"/f)

(m::Medium)(x) = refractive_index(m, x)

# * Dispersion

using Calculus

"""
    dispersion(m, d, f[, f₀=0u"Hz])

Calculate dispersion through a medium `m` of length `d`. Optionally,
remove central frequency k-vector, to keep pulse temporally centred.
"""
function dispersion(m::Medium, d::Length, f::AbstractVector{F}, f₀::Frequency = 0u"Hz") where {F<:Frequency}
    n = real(m.(f))

    # Wavevector in vacuum
    k₀ = 2π*f./u"c"
    # Wavevector in the medium
    k = n.*k₀

    if !iszero(f₀)
        # Slope of dispersion relation at central frequency
        ∂k∂f = (2π/u"c")*derivative(f -> real(m(f*u"Hz"))*f, f₀/u"Hz" .|> NoUnits)
        k -= ∂k∂f*f
    end

    exp.(-im*k*d)
end
# * Exports

export Medium, BK7, SiO2, refractive_index, dispersion
