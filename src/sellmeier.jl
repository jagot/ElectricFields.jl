# * Media

using Unitful
import Unitful: Length, Area, Frequency
import Base: show
@derived_dimension InverseArea Unitful.𝐋^(-2)

@doc raw"""
    Sellmeier(A, B, p, C, D, q)

The [Sellmeier
equations](https://en.wikipedia.org/wiki/Sellmeier_equation) are used
to describe dispersion in glass using a series of resonances, here
generalized to support variations found in the literature:

```math
n^2(\lambda) =
1 + A + \sum_i \frac {B_i\lambda^{p_i}}{\lambda^2-C_i} + \sum_j D_j \lambda^{q_j}.
```
"""
struct Sellmeier{T<:Real,Bt,Ct<:Area{T},Dt}
    A::T
    B::Vector{Bt}
    p::Vector{Int}
    C::Vector{Ct}
    D::Dt
    q::Vector{Int}
end

@doc raw"""
    Sellmeier(B, C)

Convenience constructor for [`Sellmeier`](@ref) for the canonical form of
the Sellmeier equation:

```math
n^2(\lambda) =
1 + \sum_i \frac {B_i\lambda^2}{\lambda^2-C_i}.
```
"""
Sellmeier(B::Vector{T}, C::Vector{Ct}) where {T<:Real,Ct<:Area{T}} =
    Sellmeier(zero(T), B, fill(2, length(B)), C, zeros(T,0), zeros(Int, 0))

function show(io::IO, m::Sellmeier)
    write(io, "Sellmeier(")
    show(io, m.A)
    write(io, ", ")
    show(io, m.B)
    write(io, ", ")
    show(io, m.p)
    printfmt(io, ", [{1:s}]",
             join(string.(m.C), ", "))
    printfmt(io, ", [{1:s}]",
             join(string.(m.D), ", "))
    write(io, ", ")
    show(io, m.q)
    write(io, ")")
end

Base.Broadcast.broadcastable(m::Sellmeier) = Ref(m)

# ** BK7
"""
    BK7

[Borosilicate glass](https://en.wikipedia.org/wiki/Borosilicate_glass)
is commonly used in optical lenses
"""
const BK7 = Sellmeier([1.03961212, 0.231792344, 1.01046945],
                      [6.00069867e-3u"μm"^2, 2.00179144e-2u"μm"^2, 1.03560653e2u"μm"^2])

# ** SiO₂

"""
    SiO₂

[Fused silica](https://en.wikipedia.org/wiki/Fused_quartz#Optical_properties)
"""
const SiO₂ = Sellmeier([0.6961663, 0.4079426, 0.8974794],
                       [0.0684043u"μm", 0.1162414u"μm", 9.896161u"μm"].^2)

# ** Calcite

"""
    Calcite

Calcite is a uniaxial, birefringent crystal.

Numbers taken from <http://www.newlightphotonics.com/v1/calcite-properties.html>.
"""
const Calcite = SVector(Sellmeier(1.69705, [0.0192064u"μm^2"], [0], [0.01820u"μm^2"], (-0.0151624u"μm^-2",), [2]),
                        Sellmeier(1.18438, [0.0087309u"μm^2"], [0], [0.01018u"μm^2"], (-0.0024411u"μm^-2",), [2]))

# ** Quartz

"""
    Quartz

Quartz is a uniaxial, birefrigent crystal.

Numbers taken from <https://www.newlightphotonics.com/v1/quartz-properties.html>.
"""

const Quartz = SVector(Sellmeier(1.3573, Float64[], Int[], fill(0.0u"μm^2", 0),
                                 (-0.01170u"μm^-2", 0.01054u"μm^2", 1.3414e-4u"μm^4", -4.4537e-7u"μm^6", 5.9236e-8u"μm^8"), [2,-2,-4,-6,-8] ),
                       Sellmeier(1.3849, Float64[], Int[], fill(0.0u"μm^2", 0),
                                 (-0.01259u"μm^-2", 0.01079u"μm^2", 1.6518e-4u"μm^4", -1.9474e-6u"μm^6", 9.3648e-8u"μm^8"), [2,-2,-4,-6,-8]))

# ** KTP

"""
    KTP

KTP is a biaxial crystal.

Numbers taken from <https://www.redoptronics.com/KTP-crystal.html>.
"""
const KTP = SVector(Sellmeier(1.10468, [0.89342], [2], [0.04438u"μm^2"], (-0.01036u"μm^-2",), [2]),
                    Sellmeier(1.14559, [0.87629], [2], [0.0485u"μm^2"], (-0.01173u"μm^-2",), [2]),
                    Sellmeier(0.9446,  [1.3617],  [2], [0.047u"μm^2"], (-0.01491u"μm^-2",), [2]))

# * Refractive index

function n²(m::Sellmeier, λ::Length)
    v = 1.0 + m.A
    isinf(λ) && return v

    v += NoUnits(sum([Bᵢ*λ^pᵢ/(λ^2-Cᵢ)
                     for (Bᵢ,pᵢ,Cᵢ) ∈ zip(m.B, m.p, m.C)], init=0.0))
    v += NoUnits(sum([Dⱼ*λ^qⱼ for (Dⱼ,qⱼ) ∈ zip(m.D,m.q)], init=0.0))

    v
end

n²(m::AbstractVector{<:Sellmeier}, λ::Length) = n².(m, λ)

maybe_complex(x) = x<0 ? complex(x) : x
refractive_index(m::Sellmeier, λ::Length) = √(maybe_complex(n²(m, λ)))
refractive_index(m::Sellmeier, f::Frequency) = refractive_index(m, u"c"/f)

(m::Sellmeier)(x) = refractive_index(m, x)

# * Dispersion

"""
    dispersion(m, d, f[, f₀=0u"Hz])

Calculate dispersion through a Sellmeier `m` of length `d`. Optionally,
remove central frequency k-vector, to keep pulse temporally centred.
"""
function dispersion(m::Sellmeier, d::Length, f::AbstractVector{F}, f₀::Frequency = 0u"Hz") where {F<:Frequency}
    n = real(m.(f))

    # Wavevector in vacuum
    k₀ = 2π*f./u"c"
    # Wavevector in the Sellmeier
    k = n.*k₀

    if !iszero(f₀)
        # Slope of dispersion relation at central frequency
        ∂k∂f = (2π/u"c")*ForwardDiff.derivative(f -> real(m(f*u"Hz"))*f, f₀/u"Hz" .|> NoUnits)
        k -= ∂k∂f*f
    end

    exp.(-im*k*d)
end
# * Exports

export Sellmeier, BK7, SiO₂, Calcite, Quartz, KTP, refractive_index, dispersion
