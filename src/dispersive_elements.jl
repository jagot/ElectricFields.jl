# * (An)isotropy

abstract type Tropy end
struct Isotropic <: Tropy end
struct Anisotropic <: Tropy end

Base.:(*)(::Isotropic, ::Isotropic) = Isotropic()
Base.:(*)(::Tropy, ::Tropy) = Anisotropic()

# * Dispersive elements

"""
    DispersiveElement

Base type for all dispersive elements
"""
abstract type DispersiveElement end
tropy(::DispersiveElement) = Anisotropic()
Base.Broadcast.broadcastable(de::DispersiveElement) = Ref(de)

frequency_response(::Isotropic, de::DispersiveElement, ω::AbstractVector) =
    frequency_response.(de, ω)

function frequency_response(::Anisotropic, de::DispersiveElement, ω::AbstractVector)
    T = eltype(frequency_response(de, first(ω)))
    H = zeros(T, length(ω), 3)
    for (i,ω) in enumerate(ω)
        H[i,:] .= frequency_response(de, ω)
    end
    H
end

frequency_response(de::DispersiveElement, ω::AbstractVector) =
    frequency_response(tropy(de), de, ω)

# ** Phase shift

@doc raw"""
    PhaseShift(ϕ)

Represents a phase shift according to
```math
H(\omega) =
\exp(-\im\phi).
```

# Example

```julia
julia> PhaseShift(6)
PhaseShift(ϕ = 6.0000 rad)
```
"""
struct PhaseShift{T} <: DispersiveElement
    ϕ::T
end

Base.show(io::IO, e::PhaseShift) =
    printfmt(io, "PhaseShift(ϕ = {1:0.4g} rad)", e.ϕ)

tropy(::PhaseShift) = Isotropic()

frequency_response(e::PhaseShift, ::Number) =
    exp(-im*e.ϕ)

"""
    phase_shift(f, ϕ)

Returns the field resulting from applying a [`PhaseShift`](@ref) to
the field `f`.
"""
phase_shift(f::AbstractField, ϕ; kwargs...) =
    DispersedField(f, PhaseShift(ϕ); kwargs...)

# ** Chirp

@doc raw"""
    Chirp(b, ω₀)

Represents chirp according to
```math
H(\omega) =
\exp[-\im b(ω-ω_0)^2].
```

# Example

```julia
julia> Chirp(austrip(5u"fs^2"), austrip(1.5u"eV"))
Chirp(b = 8545.5457 = 5.0000 fs², ω₀ = 0.0551 = 1.5000 eV)
```
"""
struct Chirp{T,U} <: DispersiveElement
    b::T
    ω₀::U
end

function Base.show(io::IO, c::Chirp)
    bu = u"fs^2"
    ω₀u = u"eV"
    printfmt(io, "Chirp(b = {1:0.4g} = {2:0.4g} {3:s}, ω₀ = {4:0.4g} = {5:0.4g} {6:s})",
             c.b, ustrip(auconvert(bu, c.b)), bu,
             c.ω₀, ustrip(auconvert(ω₀u, c.ω₀)), ω₀u)
end

tropy(::Chirp) = Isotropic()

frequency_response(e::Chirp, ω::Number) =
    exp(-im*e.b*(ω-e.ω₀)^2)

"""
    chirp(f, b, ω₀=photon_energy(f))

Returns the field resulting from applying a [`Chirp`](@ref) to the
field `f`.

# Example

```julia
julia> @field(F) do
           I₀ = 1.0
           T = 2.0u"fs"
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
    - A₀ = 13.1594 au
  – a Fixed carrier @ λ = 599.5849 nm (T = 2.0000 fs, ω = 0.0760 Ha = 2.0678 eV, f = 500.0000 THz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±82.68σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 58518.2144 Bohr = 3.0967 μm
  – Uₚ = 43.2922 Ha = 1.1780 keV => α = 173.1690 Bohr = 9.1637 nm

julia> chirp(F, austrip(1u"fs^2"))
DispersedField:
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
    - A₀ = 13.1594 au
  – a Fixed carrier @ λ = 599.5849 nm (T = 2.0000 fs, ω = 0.0760 Ha = 2.0678 eV, f = 500.0000 THz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±82.68σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 58518.2144 Bohr = 3.0967 μm
  – Uₚ = 43.2922 Ha = 1.1780 keV => α = 173.1690 Bohr = 9.1637 nm
  – dispersed through Chirp(b = 1709.1091 = 1.0000 fs², ω₀ = 0.0760 = 2.0678 eV)
```
"""
chirp(f::AbstractField, b, ω₀=photon_energy(f); kwargs...) =
    DispersedField(f, Chirp(b, ω₀); kwargs...)

# ** Cascade of dispersive elements

"""
    CascadedDispersiveElement(elements)

Represents the combination of multiple [`DispersiveElement`](@ref)s.

# Example

```julia
julia> PhaseShift(6)*Chirp(austrip(5u"fs^2"), austrip(1.5u"eV"))
[PhaseShift(ϕ = 6.0000 rad) ∗ Chirp(b = 8545.5457 = 5.0000 fs², ω₀ = 0.0551 = 1.5000 eV)]
```
"""
struct CascadedDispersiveElement{Elements} <: DispersiveElement
    elements::Elements
end

function Base.show(io::IO, ce::CascadedDispersiveElement)
    if isempty(ce.elements)
        write(io, "CascadedDispersiveElement(())")
    end
    write(io, "[")
    n = length(ce.elements)
    for (i,e) in enumerate(ce.elements)
        show(io, e)
        i == n || write(io, " ∗ ")
    end
    write(io, "]")
end

tropy(ce::CascadedDispersiveElement) =
    foldr((a,t) -> tropy(a) * t, ce.elements, init=Isotropic())

function frequency_response(ce::CascadedDispersiveElement, ω::Number)
    v = 1
    # We "have to" loop over the list in reverse, since the
    # multiplication rules below ensure the right-most factor in a
    # multiplication ends up as the last element of the list. Since
    # all frequency responses are scalars, they commute however.
    for e in reverse(ce.elements)
        v *= frequency_response(e, ω)
    end
    v
end

Base.:(*)(a::DispersiveElement, b::DispersiveElement) =
    CascadedDispersiveElement((a,b))

Base.:(*)(a::CascadedDispersiveElement, b::DispersiveElement) =
    CascadedDispersiveElement((a.elements...,b))

Base.:(*)(a::DispersiveElement, b::CascadedDispersiveElement) =
    CascadedDispersiveElement((a,b.elements...))

Base.:(*)(a::CascadedDispersiveElement, b::CascadedDispersiveElement) =
    CascadedDispersiveElement((a.elements...,b.elements...))

# * Find time span

function join_ranges(r1::StepRangeLen, r2::StepRangeLen)
    r1.ref == r2.ref || throw(ArgumentError("Ranges do not have the same origin"))
    step(r1) == step(r2) || throw(ArgumentError("Ranges have different steps"))

    n1 = length(r1)
    n2 = length(r2)

    StepRangeLen(r2.ref, r2.step, n1+n2, n1)
end

"""
find_overlapping_range(r1::StepRangeLen, r2::StepRangeLen)

Find the index range of `r2` which overlap with all of `r1`.
"""
function find_overlapping_range(r1::StepRangeLen, r2::StepRangeLen)
    r1.ref == r2.ref || throw(ArgumentError("Ranges do not have the same origin"))
    step(r1) == step(r2) || throw(ArgumentError("Ranges have different steps"))

    n1 = length(r1)
    n2 = length(r2)

    sel = (1:n1) .+ (r2.offset-r1.offset)
    @assert r1 == r2[sel]
    sel
end

"""
    find_time_span(f, de[, fs]; max_iter=10, ξ=2.0, tol)

Find a suitable time span ``[a,b]`` such that when `f` is dispersed
through the [`DispersiveElement`](@ref) `de`, the compact support of
the resulting field is contained in the time interval. This is done by
successively multiplying the time span on which `F` is evaluated `ξ`
before the `RFFT`, until the `IRFFT` has converged.
"""
function find_time_span(f, de, args...; max_iter=7, ξ = 2.0, tol=5e-4, verbosity=0)
    verbosity > 0 && @info "Finding large enough time span to encompass dispersed field" f de max_iter ξ tol
    a,b = endpoints(span(f))
    a == b && throw(ArgumentError("Cannot successively double an infinitesimal interval $(span(f))"))
    δt = step(timeaxis(f, args...))

    Δ = (b-a)/2
    t₀ = (a+b)/2

    tf = l -> join_ranges(reverse(range(0, stop=-l, step=-δt)), range(0, stop=l, step=δt)) .+ t₀

    t = nothing
    F = nothing

    R = Inf
    for i = 0:max_iter
        t′ = tf(Δ*ξ^i)
        if verbosity > 1 && i > 0
            println(repeat("-", 100))
            println(i)
            @show t′
        end
        H′ = frequency_response(de, rfftω(t′))
        F′ = irfft(H′.*rfft(f, t′), t′)

        if i > 0
            sel = find_overlapping_range(t, t′)
            R = norm(F - selectdim(F′, 1, sel))
            verbosity > 1 && @show R
            R < tol && break
        end

        t = t′
        F = F′
    end
    R > tol && @warn "Could not find large enough time span in $(max_iter) iterations"
    t[1] .. t[end]
end

# * Exports
export PhaseShift, Chirp, chirp
