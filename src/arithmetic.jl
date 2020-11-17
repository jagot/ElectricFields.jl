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

mutable struct SumField{A,B} <: AbstractField
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

vector_potential(f::SumField, args...) =
    vector_potential(f.a, args...) + vector_potential(f.b, args...)

polarization(f::SumField) = polarization(f.a)

function span(f::SumField)
    sa = span(f.a)
    sb = span(f.b)
    (min(sa[1], sb[1]), max(sa[2],sb[2]))
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

# ** Wrapped fields

"""
         WrappedField(field, a, b)

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

mutable struct NegatedField{F<:AbstractField} <: WrappedField
    a::F
end

(f::NegatedField)(t::Unitful.Time) = -f.a(t)
(f::NegatedField)(fs::Unitful.Frequency=default_sampling_frequency(f)) = -f.a(fs)

-(a::AbstractField,
  b::AbstractField) = a + NegatedField(b)
-(a::AbstractField) = NegatedField(a)

mutable struct NegatedCarrier <: AbstractCarrier
    carrier::AbstractCarrier
end
(carrier::NegatedCarrier)(t::Unitful.Time) = -carrier.carrier(t)

carrier(f::NegatedField) = NegatedCarrier(carrier(f.a), f.t₀)
envelope(f::NegatedField) = envelope(f.a)

Base.parent(f::NegatedField) = f.a

# ** Delayed fields

mutable struct DelayedField{F<:AbstractField,T} <: WrappedField
    a::F
    t₀::T
end

vector_potential(f::DelayedField, t) = vector_potential(f.a, t-f.t₀)

function show(io::IO, f::DelayedField)
    show(io, f.a)
    printfmt(io, "\n  – delayed by {1:.4f} jiffies = {2:s}",
             f.t₀, au2si_round(f.t₀, u"s"))
end

span(f::DelayedField) = span(parent(f)) .+ f.t₀

Base.parent(f::DelayedField) = f.a

# *** DONE Delay operators
#     Convention for delayed fields: a field delayed by a /positive/
#     time, comes /later/, i.e. we write \(f(t-\delta t)\).

delay(a::AbstractField, t₀::Unitful.Time) = DelayedField(a, austrip(t₀))
delay(a::AbstractField, nT::Real) = delay(a, nT*period(a))
delay(a::AbstractField, ϕ::Quantity{Float64, Unitful.Dimensions{()}}) = delay(a, ϕ/(2π*u"rad"))

delay(a::DelayedField) = a.t₀
delay(a::AbstractField) = 0

export delay

# ** Padded fields

"""
      PaddedField(field, a, b)

Wrapper around any electric `field`, padded with `a` units of time
before, and `b` units of time after the ordinary [`span`](@ref) of the
field.
"""
struct PaddedField{Field<:AbstractField,T} <: WrappedField
    field::Field
    a::T
    b::T

    PaddedField(field::Field, a::T, b::T) where {Field,T<:Real} =
        new{Field,T}(field, a, b)

    PaddedField(field, a::Unitful.Time, b::Unitful.Time) =
        PaddedField(field, austrip(a), austrip(b))
end

Base.parent(f::PaddedField) = f.field

function span(f::PaddedField)
    a,b = span(f.field)
    a-f.a, b+f.b
end

export PaddedField

# ** Windowed fields

"""
    WindowedField(field, a, b)

Wrapper around any electric `field`, windowed such that it is zero
outside the time interval `a..b`.
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

function span(f::WindowedField)
    a,b = span(f.field)
    max(a,f.a), min(b,f.b)
end

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
