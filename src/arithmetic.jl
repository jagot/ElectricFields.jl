# * Field arithmetic [1/1]
#   For calculations, a time-base is necessary ("normalized time"), with
#   respect to which all harmonic motions &c are analyzed.  When
#   combining fields of commensurate frequencies, it is easy to
#   establish a common time-base, most likely the fundamental wavelength
#   will be the obvious choice -- such as when adding an IR field and its
#   harmonic components as generated through e.g. HHG:
#
#   \[E(t) = \sum_q E_{2q + 1}(t)\sin[(2q+1)\omega t].\]
#
#   However, when adding incommensurate frequencies, there is no obvious
#   choice, so the user has to specify the time-base explicitly -- maybe
#   by the order in which the fields are added?
#
#   There should also be some helper routines that, when discretizing
#   the time axis, estimate whether all harmonic components of interest
#   are satisfactorily resolved.
#
# ** Sum fields

import Base: +

"""
    SumField(a, b)

The linear combination of two fields `a` and `b`.

# Example

```jldoctest
julia> @field(A) do
           I₀ = 1.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm

julia> @field(B) do
           I₀ = 0.5
           T = 1.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 5.0000e-01 au = 1.7547226e16 W cm^-2 =>
    - E₀ = 7.0711e-01 au = 363.6089 GV m^-1
    - A₀ = 0.1125 au
  – a Fixed carrier @ λ = 7.2516 nm (T = 24.1888 as, ω = 6.2832 Ha = 170.9742 eV, f = 41.3414 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±1.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 8.5598 Bohr = 452.9627 pm
  – Uₚ = 0.0032 Ha = 86.1591 meV => α = 0.0179 Bohr = 947.8211 fm

julia> A+B
┌ Linearly polarized field with
│   - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
│     - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
│     - A₀ = 0.3183 au
│   – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
│   – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
│   – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
│   – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
⊕
│ Linearly polarized field with
│   - I₀ = 5.0000e-01 au = 1.7547226e16 W cm^-2 =>
│     - E₀ = 7.0711e-01 au = 363.6089 GV m^-1
│     - A₀ = 0.1125 au
│   – a Fixed carrier @ λ = 7.2516 nm (T = 24.1888 as, ω = 6.2832 Ha = 170.9742 eV, f = 41.3414 PHz)
│   – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±1.00σ)
│   – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 8.5598 Bohr = 452.9627 pm
└   – Uₚ = 0.0032 Ha = 86.1591 meV => α = 0.0179 Bohr = 947.8211 fm
```
"""
struct SumField{A,B} <: AbstractField
    a::A
    b::B
    function SumField(a::A, b::B) where {A<:AbstractField, B<:AbstractField}
        da = dimensions(a)
        db = dimensions(b)
        da == db ||
            throw(ArgumentError("Cannot add fields of different dimensionality, got $(da) and $(db)"))
        new{A,B}(a, b)
    end
end

function show(io::IO, f::SumField)
    a_str = split(string(f.a), "\n")
    b_str = split(string(f.b), "\n")

    for (s,l) in zip("┌" * repeat("│", length(a_str)-1), a_str)
        write(io, "$s $l\n")
    end

    write(io, "⊕\n")

    for (s,l) in zip(repeat("│", length(b_str)-1) * "└", b_str)
        write(io, "$s $l\n")
    end
end

function +(a::AbstractField, b::AbstractField)
    polarization(a) == polarization(b) ||
        throw(ArgumentError("Cannot add fields of different polarization"))
    SumField(a, b)
end

for fun in [:vector_potential, :field_amplitude, :vector_potential_spectrum]
    @eval $fun(f::SumField, t::Number) =
        $fun(f.a, t) + $fun(f.b, t)
end

polarization(f::SumField) = polarization(f.a)

function span(f::SumField)
    sa = span(f.a)
    sb = span(f.b)
    min(sa.left,sb.left)..max(sa.right,sb.right)
end
steps(f::SumField, ndt::Int) = steps(f, ndt/min(austrip.(period.((f.a,f.b)))...))

for fun in [:wavelength, :period, :frequency, :wavenumber, :fundamental, :photon_energy]
    @eval begin
        function ($fun)(f::SumField)
            a = ($fun)(f.a)
            b = ($fun)(f.b)
            a != b && error("$(titlecase(string($fun))) differs between SumField composants!")
            a
        end
    end
end

max_frequency(f::SumField) =
    max(max_frequency(f.a), max_frequency(f.b))

continuity(f::SumField) =
    min(continuity(f.a), continuity(f.b))

dimensions(f::SumField) = dimensions(f.a)

phase_shift(f::SumField, ϕ) =
    SumField(phase_shift(f.a, ϕ), phase_shift(f.b, ϕ))

# ** Wrapped fields

"""
    WrappedField

Wrapper around any electric `field`
"""
abstract type WrappedField <: AbstractField end

for fun in [:params, :carrier, :envelope, :polarization,
            :wavelength, :period, :frequency, :max_frequency,
            :wavenumber, :fundamental, :photon_energy,
            :intensity, :amplitude, :duration, :continuity,
            :span, :dimensions]
    @eval $fun(f::WrappedField, args...) = $fun(parent(f), args...)
end
# [:vector_potential, :field_amplitude], should these be explicitly forwarded?

# ** Negated fields

import Base: -

"""
    NegatedField(a)

Represents a field whose [`vector_potential`](@ref) is the negative of
that of `a`.

# Example

```jldoctest
julia> @field(A) do
           I₀ = 1.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm

julia> B = -A
ElectricFields.NegatedField{ElectricFields.LinearField{ElectricFields.FixedCarrier{Quantity{Float64, 𝐋, Unitful.FreeUnits{(Eₕ^-1, ħ, c), 𝐋, nothing}}, Quantity{Float64, 𝐓, Unitful.FreeUnits{(Eₕ^-1, ħ), 𝐓, nothing}}, Float64, Int64}, ElectricFields.GaussianEnvelope{Float64}, Float64}}(Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm)

julia> field_amplitude(A, 0.5)
0.008830294133325867

julia> field_amplitude(B, 0.5)
-0.008830294133325867

julia> field_amplitude(A-A, 0.5)
0.0
```
"""
struct NegatedField{F<:AbstractField} <: WrappedField
    a::F
end

-(a::AbstractField,
  b::AbstractField) = a + NegatedField(b)
-(a::AbstractField) = NegatedField(a)

Base.parent(f::NegatedField) = f.a

for fun in [:vector_potential, :vector_potential_spectrum]
    @eval $fun(f::NegatedField, t::Number) = -$fun(parent(f), t)
end

# ** Delayed fields

"""
    DelayedField(a, t₀)

Represents a delayed copy of `a`, that appears at time `t₀` instead of
at `0`, i.e. `t₀>0` implies the field comes _later_.

# Examples

```jldoctest
julia> @field(A) do
           I₀ = 1.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm

julia> delay(A, 1u"fs")
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
  – delayed by 41.3414 jiffies = 1.0000 fs

julia> delay(A, 1.0)
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
  – delayed by 1.0000 jiffies = 24.1888 as

julia> delay(A, π*u"rad")
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm
  – delayed by 1.0000 jiffies = 24.1888 as
```
"""
struct DelayedField{F<:AbstractField,T} <: WrappedField
    a::F
    t₀::T
end

for fun in [:vector_potential, :field_amplitude, :intensity]
    @eval $fun(f::DelayedField, t::Number) =
        $fun(f.a, t-f.t₀)
end

vector_potential_spectrum(f::DelayedField, ω) =
    exp(-im*f.t₀*ω) * vector_potential_spectrum(parent(f), ω)

function show(io::IO, f::DelayedField)
    show(io, f.a)
    printfmt(io, "\n  – delayed by {1:.4f} jiffies = {2:s}",
             f.t₀, au2si_round(f.t₀, u"s"))
end

function span(f::DelayedField)
    s = span(parent(f))
    (s.left+f.t₀)..(s.right+f.t₀)
end

Base.parent(f::DelayedField) = f.a

# *** DONE Delay operators
#     Convention for delayed fields: a field delayed by a /positive/
#     time, comes /later/, i.e. we write \(f(t-\delta t)\).

delay(a::AbstractField, t₀::Union{<:Real,<:Unitful.Time}) = DelayedField(a, austrip(t₀))
delay(a::AbstractField, ϕ::Quantity{<:Real, NoDims}) = delay(a, (ϕ/(2π*u"rad"))*period(a))

delay(a::DelayedField) = a.t₀
delay(a::AbstractField) = 0

export delay

# ** Padded fields

"""
      PaddedField(field, a, b)

Wrapper around any electric `field`, padded with `a` units of time
before, and `b` units of time after the ordinary [`span`](@ref) of the
field.

# Example

```jldoctest
julia> @field(A) do
           I₀ = 1.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm

julia> B = PaddedField(A, 10.0, 30.0)
Padding before 10.0000 jiffies = 241.8884 as and after 30.0000 jiffies = 725.6653 as of
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm

julia> span(A), span(B)
(-6.0..6.0, -16.0..36.0)
```
"""
struct PaddedField{Field<:AbstractField,T} <: WrappedField
    field::Field
    a::T
    b::T

    function PaddedField(field::Field, a::T, b::T) where {Field,T<:Real}
        (a < 0 || b < 0) &&
            throw(ArgumentError("Padding must be non-negative"))
        new{Field,T}(field, a, b)
    end

    PaddedField(field, a::Unitful.Time, b::Unitful.Time) =
        PaddedField(field, austrip(a), austrip(b))
end

Base.parent(f::PaddedField) = f.field

function show(io::IO, f::PaddedField)
    printfmtln(io, "Padding before {1:.4f} jiffies = {2:s} and after {3:.4f} jiffies = {4:s} of",
               f.a, au2si_round(f.a, u"s"),
               f.b, au2si_round(f.b, u"s"))
    show(io, f.field)
end

for fun in [:vector_potential, :field_amplitude]
    @eval function $fun(f::PaddedField, t::Number)
        s = span(f.field)
        v = $fun(parent(f), clamp(t, endpoints(s)...))
        t ∈ s ? v : zero(v)
    end
end

function span(f::PaddedField)
    a,b = endpoints(span(f.field))
    (a-f.a)..(b+f.b)
end

time_integral(f::PaddedField) =
    time_integral(parent(f))

export PaddedField

# ** Windowed fields

"""
    WindowedField(field, a, b)

Wrapper around any electric `field`, windowed such that it is zero
outside the time interval `a..b`.

# Example

```jldoctest
julia> @field(A) do
           I₀ = 1.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm

julia> B = WindowedField(A, -3, 5)
Window from -3.0000 jiffies = -72.5665 as to 5.0000 jiffies = 120.9442 as of
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.3183 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm

julia> span(A), span(B)
(-6.0..6.0, -3.0..5.0)

julia> field_amplitude(A, -4)
-0.6395632362683398

julia> field_amplitude(B, -4)
0.0
```
"""
struct WindowedField{Field<:AbstractField,T} <: WrappedField
    field::Field
    a::T
    b::T

    WindowedField(field::Field, a::T, b::T) where {Field,T<:Real} =
        new{Field,T}(field, a, b)

    WindowedField(field, a::Unitful.Time, b::Unitful.Time) =
        WindowedField(field, austrip(a), austrip(b))
end

function show(io::IO, f::WindowedField)
    printfmtln(io, "Window from {1:.4f} jiffies = {2:s} to {3:.4f} jiffies = {4:s} of",
               f.a, au2si_round(f.a, u"s"),
               f.b, au2si_round(f.b, u"s"))
    show(io, f.field)
end

Base.parent(f::WindowedField) = f.field

span(f::WindowedField) = span(parent(f)) ∩ (f.a..f.b)

phase_shift(f::WindowedField, δϕ) =
    WindowedField(phase_shift(parent(f), δϕ), f.a, f.b)

for fun in [:vector_potential, :field_amplitude, :intensity]
    @eval function $fun(f::WindowedField{T}, t) where T
        v = $fun(parent(f), t)
        t < f.a || t > f.b ? zero(v) : v
    end
end

function field_amplitude(f::WindowedField, a, b)
    v = field_amplitude(parent(f), max(a, f.a), min(b, f.b))
    if a ∈ f.a..f.b || b ∈ f.a..f.b
        v
    else
        zero(v)
    end
end

function timeaxis(f::WindowedField, fs)
    t = timeaxis(f.field, fs)
    t[findall(in(f.a..f.b), t)]
end

export WindowedField
