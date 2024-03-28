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

const AnisotropicMaterial = SVector{<:Any,<:Sellmeier}

# ** BK7
@doc raw"""
    BK7

[Borosilicate glass](https://en.wikipedia.org/wiki/Borosilicate_glass)
is commonly used in optical lenses

```math
n^2(\lambda) =
1 +
\frac{1.03961212\lambda^2}{\lambda^2 - 6.00069867\times10^{-3}\;\textrm{μm}^2} +
\frac{0.231792344\lambda^2}{\lambda^2 - 2.00179144\times10^{-2}\;\textrm{μm}^2} +
\frac{1.01046945}{\lambda^2 - 1.03560653\times10^{2}\;\textrm{μm}^2}.
```
"""
const BK7 = Sellmeier([1.03961212, 0.231792344, 1.01046945],
                      [6.00069867e-3u"μm"^2, 2.00179144e-2u"μm"^2, 1.03560653e2u"μm"^2])

# ** SiO₂

@doc raw"""
    SiO₂

[Fused silica](https://en.wikipedia.org/wiki/Fused_quartz#Optical_properties)

```math
n^2(\lambda) =
1 +
\frac{0.6961663\lambda^2}{\lambda^2-(0.0684043\;\textrm{μm})^2} +
\frac{0.4079426\lambda^2}{\lambda^2-(0.1162414\;\textrm{μm})^2} +
\frac{0.8974794\lambda^2}{\lambda^2-(9.896161\;\textrm{μm})^2}.
```
"""
const SiO₂ = Sellmeier([0.6961663, 0.4079426, 0.8974794],
                       [0.0684043u"μm", 0.1162414u"μm", 9.896161u"μm"].^2)

# ** Calcite

@doc raw"""
    Calcite

Calcite is a uniaxial, birefringent crystal.

Numbers taken from <http://www.newlightphotonics.com/v1/calcite-properties.html>.

```math
\begin{aligned}
n_{\mathrm{o}}^2(\lambda)
&=
2.69705 +
\frac{0.0192064}{\lambda^2-0.01820\;\textrm{μm}^2} -
0.0151624\lambda^2, \\
n_{\mathrm{e}}^2(\lambda)
&=
2.18438 +
\frac{0.0087309}{\lambda^2-0.01018\;\textrm{μm}^2} -
0.0024411\lambda^2.
\end{aligned}
```
"""
const Calcite = SVector(Sellmeier(1.69705, [0.0192064u"μm^2"], [0], [0.01820u"μm^2"], (-0.0151624u"μm^-2",), [2]),
                        Sellmeier(1.18438, [0.0087309u"μm^2"], [0], [0.01018u"μm^2"], (-0.0024411u"μm^-2",), [2]))

# ** Quartz

@doc raw"""
    Quartz

Quartz is a uniaxial, birefrigent crystal.

Numbers taken from <https://www.newlightphotonics.com/v1/quartz-properties.html>.

```math
\begin{aligned}
n_{\mathrm{o}}^2(\lambda)
&=2.3573 -
0.01170\;\textrm{μm}^{-2}\lambda^2 +
\frac{0.01054\;\textrm{μm}^{2}}{\lambda^2} +
\frac{1.3414\times10^{-4}\;\textrm{μm}^{4}}{\lambda^4} -
\frac{4.4537\times10^{-7}\;\textrm{μm}^{6}}{\lambda^6} +
\frac{5.9236\times10^{-8}\;\textrm{μm}^{8}}{\lambda^8}, \\
n_{\mathrm{e}}^2(\lambda)
&=2.3849 -
0.01259\;\textrm{μm}^{-2}\lambda^2 +
\frac{0.01079\;\textrm{μm}^{2}}{\lambda^2} +
\frac{1.6518\times10^{-4}\;\textrm{μm}^{4}}{\lambda^4} -
\frac{1.9474\times10^{-6}\;\textrm{μm}^{6}}{\lambda^6} +
\frac{9.3648\times10^{-8}\;\textrm{μm}^{8}}{\lambda^8}.
\end{aligned}
```
"""
const Quartz = SVector(Sellmeier(1.3573, Float64[], Int[], fill(0.0u"μm^2", 0),
                                 (-0.01170u"μm^-2", 0.01054u"μm^2", 1.3414e-4u"μm^4", -4.4537e-7u"μm^6", 5.9236e-8u"μm^8"), [2,-2,-4,-6,-8] ),
                       Sellmeier(1.3849, Float64[], Int[], fill(0.0u"μm^2", 0),
                                 (-0.01259u"μm^-2", 0.01079u"μm^2", 1.6518e-4u"μm^4", -1.9474e-6u"μm^6", 9.3648e-8u"μm^8"), [2,-2,-4,-6,-8]))

# ** KTP

# """
#     KTP

# KTP is a biaxial crystal.

# Numbers taken from <https://www.newlightphotonics.com/v1/ktp-properties.html>.
# """
# const KTP = SVector(Sellmeier(2.0065, [0.03901u"μm^2"], [0], [0.04251u"μm^2"], (-0.01327u"μm^-2",), [2]),
#                     Sellmeier(2.0333, [0.04154u"μm^2"], [0], [0.04547u"μm^2"], (-0.01408u"μm^-2",), [2]),
#                     Sellmeier(2.0065, [0.05694u"μm^2"], [0], [0.05658u"μm^2"], (-0.01682u"μm^-2",), [2]))

@doc raw"""
    KTP

KTP is a biaxial crystal.

Numbers taken from <https://www.redoptronics.com/KTP-crystal.html>.

```math
\begin{aligned}
n_x^2(\lambda)
&= 2.10468 +
\frac{0.89342\lambda^2}{\lambda^2-0.04438\;\textrm{μm}^2} -
0.01036\;\textrm{μm}^{-2}\lambda^2, \\
n_y^2(\lambda)
&= 2.14559 +
\frac{0.87629\lambda^2}{\lambda^2-0.0485\;\textrm{μm}^2} -
0.01173\;\textrm{μm}^{-2}\lambda^2, \\
n_z^2(\lambda)
&= 1.9446 +
\frac{1.3617\lambda^2}{\lambda^2-0.047\;\textrm{μm}^2} -
0.01491\;\textrm{μm}^{-2}\lambda^2. \\
\end{aligned}
```
"""
const KTP = SVector(Sellmeier(1.10468, [0.89342], [2], [0.04438u"μm^2"], (-0.01036u"μm^-2",), [2]),
                    Sellmeier(1.14559, [0.87629], [2], [0.0485u"μm^2"], (-0.01173u"μm^-2",), [2]),
                    Sellmeier(0.9446,  [1.3617],  [2], [0.047u"μm^2"], (-0.01491u"μm^-2",), [2]))

# * Refractive index

function n²(m::Sellmeier, λ::Length)
    O = one(ustrip(λ))
    v = O + m.A
    isinf(λ) && return v

    v += sum([NoUnits(Bᵢ*λ^pᵢ/(λ^2-Cᵢ))
             for (Bᵢ,pᵢ,Cᵢ) ∈ zip(m.B, m.p, m.C)], init=zero(O))
    v += sum([NoUnits(Dⱼ*λ^qⱼ)
             for (Dⱼ,qⱼ) ∈ zip(m.D,m.q)], init=zero(O))

    v
end

n²(m::AbstractVector{<:Sellmeier}, λ::Length) = n².(m, λ)

maybe_complex(x) = x<0 ? complex(x) : x
refractive_index(m::Sellmeier, λ::Length) = √(maybe_complex(n²(m, λ)))
refractive_index(m::Sellmeier, f::Frequency) = refractive_index(m, u"c"/f)

(m::Sellmeier)(x) = refractive_index(m, x)

# * Dispersion

function dispersion_slope(m::Sellmeier, ω₀)
    f₀ = ustrip(auconvert(u"Hz", ω₀/2π))
    inv(austrip(1u"c"))*ForwardDiff.derivative(f -> real(m(f*u"Hz"))*f, f₀)
end

# We can only remove a common slope, since the different axes _should_
# be delayed with respect to each other.
dispersion_slope(m::AnisotropicMaterial, ω₀) =
    mean(dispersion_slope.(m, ω₀))

# * Exports

export Sellmeier, BK7, SiO₂, Calcite, Quartz, KTP, refractive_index, dispersion, dispersion_slope
