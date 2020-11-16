# * TODO Envelopes [2/3]
#   The envelopes implemented below are all /amplitude/ envelopes,
#   since that is what is being used in calculations. However, they may
#   be specified using intensity-related quantities, e.g. Gaussian
#   pulses are most often specified using the FWHM duration of their
#   /intensity/ envelopes.
#

envelope_types = Dict{Symbol,Any}()

# ** DONE Gaussian
@doc raw"""
    GaussianEnvelope

A Gaussian pulse is given by

```math
I_0\exp\left(-\frac{t^2}{2\sigma^2}\right),
```

where the standard deviation ``σ`` is related to the FWHM duration τ
of the intensity envelope as

```math
\sigma = \frac{\tau}{2\sqrt{2\ln 2}}.
```

Furthermore, the _amplitude_ standard deviation ``σ′`` is proportional
to the intensity ditto: ``\sigma' = \sqrt{2}\sigma``. Therefore, the
amplitude envelope is given by

```math
E_0\exp\left(-\frac{t^2}{2{\sigma'}^2}\right)
=E_0\exp\left(-\frac{t^2}{4\sigma^2}\right)
=E_0\exp\left(-\frac{2\ln2t^2}{\tau^2}\right).
```

Since a Gaussian never ends, we specify how many ``σ`` we require; the
resulting time window will be rounded up to an integer amount of
cycles of the fundamental.
"""
struct GaussianEnvelope{T} <: AbstractEnvelope
    τ::T # Intensity FWHM
    α::T # Vector potential exponential coefficient
    σmax::T
    tmax::T # Maximum time. Time window: [-tmax,tmax]
end
envelope_types[:gauss] = GaussianEnvelope

(env::GaussianEnvelope)(t) = exp(-env.α*t^2)

show(io::IO, env::GaussianEnvelope) =
    printfmt(io, "Gaussian envelope of duration {1:s} (intensity FWHM; ±{2:0.2f}σ) ",
             au2si_round(env.τ, u"s"), env.σmax)

function gaussian_common!(field_params, carrier; verbosity=0)
    @namespace!(field_params) do
        if τ
            σ = τ/(2*√(2log(2)))
        else
            τ = 2*√(2log(2))*σ
        end

        if σmax
            Tmax = ceil(Int, σmax*σ/T)
            tmax = Tmax*T
        else
            if tmax
                Tmax = ceil(Int, tmax/T)
            elseif Tmax
                tmax = Tmax*T
            end
            tmax = Tmax*T
        end
        σmax = tmax/σ
    end
    τ = austrip(field_params[:τ])
    α = 2log(2)/τ^2
    period = austrip(field_params[:T])
    ω = austrip(field_params[:ω])

    σmax = austrip(field_params[:σmax])
    tmax = austrip(field_params[:tmax])

    τ < 0.5period &&
        @warn "Pulse durations smaller than half a period not reliably supported"

    # Find the exponential coefficient for the vector potential that
    # will give the desired FWHM for the intensity profile.
    f = α -> begin
        env = GaussianEnvelope(τ, α[1], σmax, tmax)
        # TODO: This is not valid for e.g. elliptical fields, since
        # the intensity is given as a sum(abs2, ...) over the z- and
        # x-components, which will conceivably require another
        # exponential coefficient for the vector potential to yield
        # the desired FWHM.
        field = make_temp_field(carrier, env, field_params)
        abs2(intensity(field, τ/2) - ω^2/2)
    end
    verbosity > 0 && @info "Finding optimal vector potential coefficient"
    res = optimize(f, [α], BFGS())
    verbosity > 0 && display(res)

    field_params[:α] = res.minimizer[1]
end

function GaussianEnvelope(field_params::Dict{Symbol,Any}, carrier)
    test_field_parameters(field_params, [:T]) # Period time required to round time window up
    test_field_parameters(field_params, [:τ, :σ])
    test_field_parameters(field_params, [:σmax, :tmax, :Tmax])

    gaussian_common!(field_params, carrier)

    @unpack τ, α, σmax, tmax = field_params
    GaussianEnvelope(austrip(τ), α, σmax, austrip(tmax))
end

continuity(::GaussianEnvelope) = Inf
span(env::GaussianEnvelope) = (-env.tmax, env.tmax)

time_integral(env::GaussianEnvelope) = env.σ*√(2π)

# *** Spectrum
@doc raw"""
    spectrum(env::GaussianEnvelope)

Gaussians belong to the [Schwartz
class](https://en.wikipedia.org/wiki/Schwartz_space), i.e. functions
who, under Fourier transform, are mapped back to the same space. That
is to say, the Fourier transform of a Gaussian is a Gaussian:

```math
\exp(-\alpha t^2) \leftrightarrow
\frac{1}{\sqrt{2\alpha}}
\exp\left(-\frac{\omega^2}{4\alpha}\right).
```

Comparing with the above, we find that the spectral standard
deviation

```math
\Omega = \sqrt{2\alpha} = \frac{2\sqrt{\ln 2}}{\tau},
```

and the Gaussian function in the spectral domain is thus

```math
E(\omega) =
\frac{E_0\tau}{2\sqrt{\ln 2}}
\exp\left[-\frac{(\omega\tau)^2}{8\ln2}\right].
```
"""
function spectrum(env::GaussianEnvelope)
    N = env.E₀*env.τ/(2*√(log(2)))
    ω -> N*exp(-(ω*env.τ)^2/8log(2))
end

# ** Truncated Gaussian

struct TruncatedGaussianEnvelope{T} <: AbstractEnvelope
    τ::T # Intensity FWHM
    σ::T # Intensity std.dev.
    α::T # Vector potential exponential coefficent
    toff::T # Start of hard turn-off
    tmax::T # Maximum time; end of hard turn-off. Time window: [-tmax,tmax]
    Tmax::Int # Maximum time, in cycles of the fundamental.
end
envelope_types[:trunc_gauss] = TruncatedGaussianEnvelope

function (env::TruncatedGaussianEnvelope{T})(t) where T
    at = abs(t)
    at > env.tmax && return zero(T)
    @unpack α, toff, tmax = env
    if at ≤ env.toff
        exp(-α*t^2)
    else
        exp(-α*(toff + 2/π*(tmax-toff)*tan(π/2*(at-toff)/(tmax-toff)))^2)
    end
end

show(io::IO, env::TruncatedGaussianEnvelope) =
    printfmt(io, "Truncated Gaussian envelope of duration {1:s} jiffies = {2:s} (intensity FWHM; turn-off from {3:s} to {4:s}) ",
             env.τ, au2si_round(env.τ, u"s"),
             au2si_round(env.toff, u"s"), au2si_round(env.tmax, u"s"))

function TruncatedGaussianEnvelope(field_params::Dict{Symbol,Any}, carrier)
    test_field_parameters(field_params, [:T]) # Period time required to round time window up
    test_field_parameters(field_params, [:τ, :σ])
    test_field_parameters(field_params, [:toff])
    test_field_parameters(field_params, [:σmax, :tmax, :Tmax])

    gaussian_common!(field_params, carrier)

    @unpack τ, σ, α, toff, tmax, Tmax = field_params
    toff < tmax || throw(ArgumentError("Hard turn-off must occur before end of pulse"))
    TruncatedGaussianEnvelope(austrip(τ), austrip(σ), α, austrip(toff), austrip(tmax), Tmax)
end

continuity(::TruncatedGaussianEnvelope) = Inf # This is not exactly true
span(env::TruncatedGaussianEnvelope) = (-env.tmax, env.tmax)

# TODO: Take truncation into account
time_integral(env::TruncatedGaussianEnvelope) = env.σ*√(2π)

# ** DONE Trapezoidal

struct TrapezoidalEnvelope{T} <: AbstractEnvelope
    ramp_up::T
    flat::T
    ramp_down::T
    period::T
end
envelope_types[:trapezoidal] = TrapezoidalEnvelope
# Common alias
envelope_types[:tophat] = TrapezoidalEnvelope

function (env::TrapezoidalEnvelope{T})(t) where T
    t /= env.period
    if t < 0
        zero(T)
    elseif t < env.ramp_up
        t/env.ramp_up
    elseif t < env.ramp_up + env.flat
        1
    elseif t < env.ramp_up + env.flat + env.ramp_down
        1 - (t-env.ramp_up-env.flat)/env.ramp_down
    else
        zero(T)
    end
end

show(io::IO, env::TrapezoidalEnvelope) =
    printfmt(io, "/{1:d}‾{2:d}‾{3:d}\\ cycles trapezoidal envelope",
             env.ramp_up, env.flat, env.ramp_down)

function TrapezoidalEnvelope(field_params::Dict{Symbol,Any}, args...)
    test_field_parameters(field_params, [:T]) # Period time required to relate ramps/flat to cycles
    test_field_parameters(field_params, [:ramp, :ramp_up])
    test_field_parameters(field_params, [:ramp, :ramp_down])
    test_field_parameters(field_params, [:flat])

    @namespace!(field_params) do
        if ramp
            ramp_up = ramp
            ramp_down = ramp
        end
    end

    @unpack ramp_up, flat, ramp_down, T = field_params

    ramp_up >= 0 || error("Negative up-ramp not supported")
    flat >= 0 || error("Negative flat region not supported")
    ramp_down >= 0 || error("Negative down-ramp not supported")
    ramp_up + flat + ramp_down > 0 || error("Pulse length must be non-zero")

    TrapezoidalEnvelope(ramp_up, flat, ramp_down, austrip(T))
end

continuity(::TrapezoidalEnvelope) = 0
span(env::TrapezoidalEnvelope{T}) where T =
    (zero(T), (env.ramp_up + env.flat + env.ramp_down)*env.period)

# ** Cos²

struct Cos²Envelope{T} <: AbstractEnvelope
    cycles::T
    T::T
end
envelope_types[:cos²] = Cos²Envelope
envelope_types[:cos2] = Cos²Envelope

function (env::Cos²Envelope{T})(t) where T
    t /= (env.cycles*env.T)
    if -1 ≤ 2t ≤ 1
        cospi(t)^2
    else
        zero(T)
    end
end

show(io::IO, env::Cos²Envelope) =
    printfmt(io, "{1:0.2f} cycles cos² envelope", env.cycles)

function Cos²Envelope(field_params::Dict{Symbol,Any}, args...)
    test_field_parameters(field_params, [:T]) # Period time required to relate ramps/flat to cycles
    test_field_parameters(field_params, [:cycles])

    @unpack cycles, T = field_params

    cycles >= 0 || error("Negative duration not supported")

    Cos²Envelope(cycles, austrip(T))
end

continuity(::Cos²Envelope) = 0
function span(env::Cos²Envelope)
    s = env.cycles*env.T/2
    (-s, s)
end

# ** Exports

export continuity
