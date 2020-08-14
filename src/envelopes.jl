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
struct GaussianEnvelope <: AbstractEnvelope
    τ::Number # Intensity FWHM
    σ::Number # Intensity std.dev.
    σ′::Number # Envelope std.dev.
    σmax::Number
    σ′max::Number
    tmax::Number # Maximum time. Time window: [-tmax,tmax]
    Tmax::Integer # Maximum time, in cycles of the fundamental.
    I₀::Number
    E₀::Number
end
envelope_types[:gauss] = GaussianEnvelope

(env::GaussianEnvelope)(t::Unitful.Time) = env.E₀*exp(-t^2/(2*env.σ′^2))

show(io::IO, env::GaussianEnvelope) =
    write(io, @sprintf("I₀ = %0.2g %s Gaussian envelope of duration %0.2g %s (intensity FWHM; ±%0.2fσ) ",
                       usplit(env.I₀)..., usplit(env.τ)..., env.σmax))

function gaussian_common!(field_params)
    @namespace!(field_params) do
        if τ
            σ = τ/(2*√(2log(2)))
        else
            if σ′
                σ = σ′/√2
            end
            τ = 2*√(2log(2))*σ
        end
        if !σ′
            σ′ = √2*σ
        end

        if σmax || σ′max
            if σmax
                Tmax = ceil(Int, σmax*σ/T)
            elseif σ′max
                Tmax = ceil(Int, σ′max*σ′/T)
            end
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
        σ′max = tmax/σ′
    end
end

function GaussianEnvelope(field_params::Dict{Symbol,Any})
    test_field_parameters(field_params, [:T]) # Period time required to round time window up
    test_field_parameters(field_params, [:τ, :σ, :σ′])
    test_field_parameters(field_params, [:σmax, :σ′max, :tmax, :Tmax])

    gaussian_common!(field_params)

    @unpack τ, σ, σ′, σmax, σ′max, tmax, Tmax, I₀, E₀ = field_params
    GaussianEnvelope(τ, σ, σ′, σmax, σ′max, tmax, Tmax, I₀, E₀)
end

continuity(::GaussianEnvelope) = Inf
span(env::GaussianEnvelope) = (-env.tmax, env.tmax)

intensity(env::GaussianEnvelope) = env.I₀
amplitude(env::GaussianEnvelope) = env.E₀

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

struct TruncatedGaussianEnvelope <: AbstractEnvelope
    τ::Unitful.Time # Intensity FWHM
    σ::Unitful.Time # Intensity std.dev.
    σ′::Unitful.Time # Envelope std.dev.
    toff::Unitful.Time # Start of hard turn-off
    tmax::Unitful.Time # Maximum time; end of hard turn-off. Time window: [-tmax,tmax]
    Tmax::Integer # Maximum time, in cycles of the fundamental.
    I₀::Intensity
    E₀::ElectricField
end
envelope_types[:trunc_gauss] = TruncatedGaussianEnvelope

function (env::TruncatedGaussianEnvelope)(t::Unitful.Time)
    at = abs(t)
    at > env.tmax && return zero(env.E₀)
    α = 1/(2env.σ′^2)
    amp = if at ≤ env.toff
        exp(-α*t^2)
    else
        exp(-α*(env.toff + 2/π*(env.tmax-env.toff)*tan(π/2*(at-env.toff)/(env.tmax-env.toff)))^2)
    end
    env.E₀*amp
end

show(io::IO, env::TruncatedGaussianEnvelope) =
    @printf(io, "I₀ = %0.2g %s truncated Gaussian envelope of duration %0.2g %s (intensity FWHM; turn-off from %0.2g %s to %0.2g %s) ",
            usplit(env.I₀)..., usplit(env.τ)..., usplit(env.toff)..., usplit(env.tmax)...)

function TruncatedGaussianEnvelope(field_params::Dict{Symbol,Any})
    test_field_parameters(field_params, [:T]) # Period time required to round time window up
    test_field_parameters(field_params, [:τ, :σ, :σ′])
    test_field_parameters(field_params, [:toff])
    test_field_parameters(field_params, [:σmax, :σ′max, :tmax, :Tmax])

    gaussian_common!(field_params)

    @unpack τ, σ, σ′, toff, tmax, Tmax, I₀, E₀ = field_params
    toff < tmax || throw(ArgumentError("Hard turn-off must occur before end of pulse"))
    TruncatedGaussianEnvelope(τ, σ, σ′, toff, tmax, Tmax, I₀, E₀)
end

continuity(::TruncatedGaussianEnvelope) = Inf # This is not exactly true
span(env::TruncatedGaussianEnvelope) = (-env.tmax, env.tmax)

intensity(env::TruncatedGaussianEnvelope) = env.I₀
amplitude(env::TruncatedGaussianEnvelope) = env.E₀

# ** DONE Trapezoidal

struct TrapezoidalEnvelope <: AbstractEnvelope
    ramp_up::Number
    flat::Number
    ramp_down::Number
    I₀::Number
    E₀::Number
    T::Unitful.Time
end
envelope_types[:trapezoidal] = TrapezoidalEnvelope
# Common alias
envelope_types[:tophat] = TrapezoidalEnvelope

function (env::TrapezoidalEnvelope)(t::Unitful.Time)
    t /= env.T
    f = if t < 0
        0
    elseif t < env.ramp_up
        t/env.ramp_up
    elseif t < env.ramp_up + env.flat
        1
    elseif t < env.ramp_up + env.flat + env.ramp_down
        1 - (t-env.ramp_up-env.flat)/env.ramp_down
    else
        0
    end
    f*env.E₀
end

show(io::IO, env::TrapezoidalEnvelope) =
    write(io, @sprintf("I₀ = %0.2g %s /%s‾%s‾%s\\ cycles trapezoidal envelope",
                       usplit(env.I₀)...,
                       env.ramp_up, env.flat, env.ramp_down))

function TrapezoidalEnvelope(field_params::Dict{Symbol,Any})
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

    @unpack ramp_up, flat, ramp_down, I₀, E₀, T = field_params

    ramp_up >= 0 || error("Negative up-ramp not supported")
    flat >= 0 || error("Negative flat region not supported")
    ramp_down >= 0 || error("Negative down-ramp not supported")
    ramp_up + flat + ramp_down > 0 || error("Pulse length must be non-zero")

    TrapezoidalEnvelope(ramp_up, flat, ramp_down, I₀, E₀, T)
end

continuity(::TrapezoidalEnvelope) = 0
span(env::TrapezoidalEnvelope) = (0u"fs", (env.ramp_up + env.flat + env.ramp_down)*env.T)

intensity(env::TrapezoidalEnvelope) = env.I₀
amplitude(env::TrapezoidalEnvelope) = env.E₀

# ** TODO Sin2
#
# ** CW

struct CWEnvelope <: AbstractEnvelope
    tmax::Number # Maximum time. Time window: [-tmax,tmax]
    Tmax::Integer # Maximum time, in cycles of the fundamental.
    I₀::Intensity
    E₀::ElectricField
end
envelope_types[:cw] = CWEnvelope

(env::CWEnvelope)(t::Unitful.Time) = env.E₀

show(io::IO, env::CWEnvelope) =
    write(io, @sprintf("I₀ = %0.2g %s CW envelope of duration %0.2g %s (%0.2g cycles) ",
                       usplit(env.I₀)..., usplit(env.tmax)..., env.Tmax))

function CWEnvelope(field_params::Dict{Symbol,Any})
    test_field_parameters(field_params, [:T]) # Period time required to set time window
    test_field_parameters(field_params, [:tmax, :Tmax])

    @namespace!(field_params) do
        if tmax
            Tmax = ceil(Int, tmax/T)
        end
        tmax = Tmax*T
    end

    @unpack tmax, Tmax, I₀, E₀ = field_params
    CWEnvelope(tmax, Tmax, I₀, E₀)
end

continuity(::CWEnvelope) = Inf
span(env::CWEnvelope) = (0u"s", env.tmax)

intensity(env::CWEnvelope) = env.I₀
amplitude(env::CWEnvelope) = env.E₀

# ** Exports

export continuity
