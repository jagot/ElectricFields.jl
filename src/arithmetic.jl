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

mutable struct SumField <: AbstractField
    a::AbstractField
    b::AbstractField
end

function show(io::IO, f::SumField)
    a_str = split(string(f.a), "\n")
    b_str = split(string(f.b), "\n")

    for (s,l) in zip("⌈" * repeat("|", length(a_str)-1), a_str)
        write(io, "$s $l\n")
    end

    write(io, "⊕\n")

    for (s,l) in zip(repeat("|", length(b_str)-1) * "⌊", b_str)
        write(io, "$s $l\n")
    end
end

+(a::AbstractField,
  b::AbstractField) = SumField(a, b)

(f::SumField)(t::Unitful.Time) = f.a(t) + f.b(t)
(f::SumField)(fs::Unitful.Frequency=default_sampling_frequency(f)) = f.a(fs) + f.b(fs)

function span(f::SumField)
    sa = span(f.a)
    sb = span(f.b)
    (min(sa[1], sb[1]), max(sa[2],sb[2]))
end

for fun in [:wavelength, :period, :frequency, :wavenumber, :fundamental, :photon_energy]
    @eval begin
        function ($fun)(f::SumField)
            a = ($fun)(f.a)
            b = ($fun)(f.b)
            a != b && error("$(ucfirst(string($fun))) differs between SumField composants!")
            a
        end
    end
end

max_frequency(f::SumField) =
    max(max_frequency(f.a), max_frequency(f.b))

continuity(f::SumField) =
    min(continuity(f.a), continuity(f.b))

# ** Negated fields

import Base: -

mutable struct NegatedField <: AbstractField
    a::AbstractField
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

# ** Delayed fields

mutable struct DelayedField <: AbstractField
    a::AbstractField
    t₀::Number
end
(f::DelayedField)(t::Unitful.Time) = f.a(t-f.t₀)

function show(io::IO, f::DelayedField)
    show(io, f.a)
    write(io, "\n  – delayed by ")
    show(io, f.t₀)
end

mutable struct DelayedCarrier <: AbstractCarrier
    carrier::AbstractCarrier
    t₀::Number
end
(carrier::DelayedCarrier)(t::Unitful.Time) = carrier.carrier(t-carrier.t₀)

mutable struct DelayedEnvelope <: AbstractEnvelope
    env::AbstractEnvelope
    t₀::Number
end
(envelope::DelayedEnvelope)(t::Unitful.Time) = envelope.env(t-envelope.t₀)

carrier(f::DelayedField) = DelayedCarrier(carrier(f.a), f.t₀)
envelope(f::DelayedField) = DelayedEnvelope(envelope(f.a), f.t₀)

span(env::DelayedEnvelope) = span(env.env) .+ env.t₀

for FieldType in [:NegatedField, :DelayedField]
    for fun in [:wavelength, :period, :frequency, :max_frequency,
                :wavenumber, :fundamental, :photon_energy,
                :intensity, :amplitude, :duration, :continuity,
                :span, :steps]
        @eval ($fun)(f::($FieldType)) = ($fun)(f.a)
    end
end

# *** DONE Delay operators
#     Convention for delayed fields: a field delayed by a /positive/
#     time, comes /later/, i.e. we write \(f(t-\delta t)\).

delay(a::AbstractField, t₀::Unitful.Time) = DelayedField(a, t₀)
delay(a::AbstractField, nT::Real) = delay(a, nT*period(a))
delay(a::AbstractField, ϕ::Quantity{Float64, Unitful.Dimensions{()}}) = delay(a, ϕ/(2π*u"rad"))

delay(a::DelayedField) = a.t₀
delay(a::AbstractField) = 0u"s"

export delay

# ** Wrapped fields

"""
         WrappedField(field, a, b)

     Wrapper around any electric `field`
"""
abstract type WrappedField <: AbstractField end

for fun in [:params, :carrier, :envelope, :span, :vector_potential, :field_amplitude, :amplitude, :intensity]
    @eval $fun(f::WrappedField, args...) = $fun(parent(f), args...)
end

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
