# * Dispersed fields
#   We calculate the effect of dispersion described by the Sellmeier
#   equations in the frequency domain, to which we transform via an
#   FFT. For this to be possible, we only allow the evaluation of the
#   field for a specified sampling frequency, i.e. we don't provide an
#   implementation for evaluating a =DispersedField= at an arbitrary
#   time point.

struct DispersedField <: AbstractField
    a::AbstractField
    m::Medium
    d::Unitful.Length
end

function show(io::IO, f::DispersedField)
    show(io, f.a)
    write(io, "\n  – dispersed through $(f.d) of $(f.m)")
end

disperse(a::AbstractField, m::Medium, d::Unitful.Length) =
    DispersedField(a, m, d)

(f::DispersedField)(t::Union{Unitful.Time,AbstractVector{<:Unitful.Time}}) =
    error("Dispersed fields can only be evaluated by specifying a sampling frequency")

function (f::DispersedField)(fs::Unitful.Frequency = default_sampling_frequency(f.a))
    F̂ = fft(f.a, fs)
    ifft(F̂.*dispersion(f.m, f.d, fftfreq(f.a, fs), frequency(f.a)))*(u"V"*u"m")
end

max_frequency(f::DispersedField) = max_frequency(f.a)

# This is, strictly speaking, not true, since dispersing a pulse will
# in general stretch it in the time domain. It is up to the user to
# ensure that the time window is large enough to contain even the
# stretched pulse, by setting the σmax (or friends) to a large enough
# value.
span(f::DispersedField) = span(f.a)

export disperse
