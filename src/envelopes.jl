# * Envelopes
#   The envelopes implemented below are all /vector potential/
#   envelopes, since that is what is being used in
#   calculations. However, they may be specified using
#   intensity-related quantities, e.g. Gaussian pulses are most often
#   specified using the FWHM duration of their /intensity/ envelopes.
#

envelope_types = Dict{Symbol,Any}()

complex_squared(t) = t*t

# ** Gaussian-like
abstract type AbstractGaussianEnvelope{T} <: AbstractEnvelope end

# *** Gaussian

@doc raw"""
    GaussianEnvelope

A Gaussian pulse is given by

```math
f(t) = \exp\left(-\frac{t^2}{2\sigma^2}\right),
```

where the standard deviation ``σ`` is related to the FWHM duration τ
of the intensity envelope as

```math
\sigma = \frac{\tau}{2\sqrt{2\ln 2}}.
```

Since we define all fields in terms of the vector potential ``A(t)``,
and the [`instantaneous_intensity`](@ref) is given by
``\abs{-\partial_t A(t)}^2``, we have to find a coefficient ``α`` such
that the field amplitude

```math
F(t) \sim -\partial_t \exp(-\alpha t^2) \sin(\omega t + \phi)
```

has an _intensity envelope_ with the desired FWHM; we do this
iteratively in [`gaussian_common!`](@ref). This is mainly important
for ultrashort pulses, since for long pulse durations we can
approximate ``\exp(-\alpha t^2) \sim 1`` such that only the carrier
contributes to the time derivative.

Since a Gaussian never ends, we specify how many ``σ`` we require; the
resulting time window will be rounded up to an integer amount of
cycles of the fundamental.

Required parameters:
- `λ|T|f|ν|ω|ħω`,
- `τ|σ|α`,
- `σmax|tmax|Tmax`,
- `env=:gauss` (optional, since [`GaussianEnvelope`](@ref) is the default envelope).
"""
struct GaussianEnvelope{T} <: AbstractGaussianEnvelope{T}
    τ::T # Intensity FWHM
    σ::T # Intensity std.dev.
    α::T # Vector potential exponential coefficient
    σmax::T
    tmax::T # Maximum time. Time window: [-tmax,tmax]
end
envelope_types[:gauss] = GaussianEnvelope

(env::GaussianEnvelope)(t) = exp(-env.α*complex_squared(t))

show(io::IO, env::GaussianEnvelope) =
    printfmt(io, "Gaussian envelope of duration {1:s} (intensity FWHM; ±{2:0.2f}σ)",
             au2si_round(env.τ, u"s"), env.σmax)

function findα(τ, ω, ϵ=0.2)
    eq(α) = 2log(2)/τ^2*(1 + 1/log(2)*log(1 + τ^2*α^2/ω^2)) - α

    α₀ = 2log(2)/τ^2
    x = τ*ω/2π
    # This is our heuristic upper bound
    w = 1 + ϵ + 4log(1 + 2/x)*exp(-x)
    find_zero(eq, (α₀,w*α₀))
end

function findτ(α, ω)
    eq(τ) = 2log(2)/τ^2*(1 + 1/log(2)*log(1 + τ^2*α^2/ω^2)) - α

    τ₀ = √(2log(2)/α)
    x = α/ω^2
    # This is our heuristic upper bound
    w̃ = 1 + 85sqrt(log(1 + √x))/100 * (2 - exp(-x))
    find_zero(eq, (τ₀,w̃*τ₀))
end

@doc raw"""
    gaussian_common!(field_params, carrier[; Tmax_rounder, verbosity])

Compute parameters common to Gaussian envelopes,
i.e. [`GaussianEnvelope`](@ref) and
[`TruncatedGaussianEnvelope`](@ref); most importantly, given an
intensity FWHM or σ, we need to figure out the coefficient ``\alpha``
for the envelope of the vector potential ``\exp(-\alpha t^2)``, such
that the electric field amplitude and intensity have the desired
durations. The optional function `Tmax_rounder` determines if the
envelop should be extended to encompass an integer amount of cycles of
the `carrier` (default).
"""
function gaussian_common!(field_params, carrier;
                          Tmax_rounder = Base.Fix1(ceil, Int), verbosity=0)
    field_params[:T] = austrip(field_params[:T])
    field_params[:ω] = austrip(field_params[:ω])
    @namespace!(field_params) do
        if α
            τ = findτ(α, ω)
        end

        if τ
            τ = austrip(τ)
            σ = τ/(2*√(2log(2)))
        else
            σ = austrip(σ)
            τ = 2*√(2log(2))*σ
        end

        if σmax
            Tmax = Tmax_rounder(σmax*σ/T)
            tmax = Tmax*T
        else
            if tmax
                tmax = austrip(tmax)
                Tmax = Tmax_rounder(tmax/T)
            elseif Tmax
                Tmax = austrip(Tmax)
                tmax = Tmax*T
            end
            tmax = Tmax*T
        end
        σmax = tmax/σ
    end
    τ = austrip(field_params[:τ])
    if :α ∉ keys(field_params)
        field_params[:α] = findα(τ, ω)
    end
end

function GaussianEnvelope(field_params::Dict{Symbol,Any}, carrier)
    test_field_parameters(field_params, [:T]) # Period time required to round time window up
    test_field_parameters(field_params, [:τ, :σ, :α])
    test_field_parameters(field_params, [:σmax, :tmax, :Tmax])

    gaussian_common!(field_params, carrier)

    @unpack τ, σ, α, σmax, tmax = field_params
    GaussianEnvelope(austrip(τ), austrip(σ), α, austrip(σmax), austrip(tmax))
end

duration(env::GaussianEnvelope) = env.τ
continuity(::GaussianEnvelope) = Inf
span(env::GaussianEnvelope) = -env.tmax..env.tmax

time_integral(env::GaussianEnvelope) = env.σ*√(2π)

# *** Truncated Gaussian

@doc raw"""
    TruncatedGaussianEnvelope

Since a Gaussian function never ends, this envelope adds a soft truncation over a time interval:

```math
f(t)=
\begin{cases}
\exp(-\alpha t^2), & \abs{t} \leq t_{\textrm{off}},\\
\exp\left\{
-\alpha\left[
t_{\textrm{off}} + \frac{2}{\pi}(t_{\textrm{max}} - t_{\textrm{off}})
\tan\left(
\frac{\pi}{2}
\frac{\abs{t} - t_{\textrm{off}}}{t_{\textrm{max}} - t_{\textrm{off}}}
\right)
\right]^2
\right\},
& t_{\textrm{off}} < \abs{t} \leq t_{\textrm{max}}, \\
0, & t_{\textrm{max}} < \abs{t}.
\end{cases}
```
This is Eq. (72) of

- Patchkovskii, S., & Muller, H. (2016). Simple, accurate, and
  efficient implementation of 1-electron atomic time-dependent
  Schrödinger equation in spherical coordinates. Computer Physics
  Communications, 199,
  153–169. [10.1016/j.cpc.2015.10.014](http://dx.doi.org/10.1016/j.cpc.2015.10.014)

Required parameters:
- `λ|T|f|ν|ω|ħω`,
- `τ|σ|α`,
- `toff`,
- `σmax|tmax|Tmax`,
`env=:trunc_gauss`.

Beyond this, everything else is the same as for
[`GaussianEnvelope`](@ref).
"""
struct TruncatedGaussianEnvelope{T} <: AbstractGaussianEnvelope{T}
    τ::T # Intensity FWHM
    σ::T # Intensity std.dev.
    α::T # Vector potential exponential coefficent
    toff::T # Start of hard turn-off
    tmax::T # Maximum time; end of hard turn-off. Time window: [-tmax,tmax]
    Tmax::T # Maximum time, in cycles of the fundamental.
end
envelope_types[:trunc_gauss] = TruncatedGaussianEnvelope

function (env::TruncatedGaussianEnvelope{T})(t) where T
    at = abs(real(t))
    at > env.tmax && return zero(T)
    @unpack α, toff, tmax = env
    if at ≤ env.toff
        exp(-α*complex_squared(t))
    else
        exp(-α*complex_squared(toff + 2/π*(tmax-toff)*tan(π/2*(at-toff)/(tmax-toff))))
    end
end

show(io::IO, env::TruncatedGaussianEnvelope) =
    printfmt(io, "Truncated Gaussian envelope of duration {1:.4f} jiffies = {2:s} (intensity FWHM; turn-off from {3:s} to {4:s})",
             env.τ, au2si_round(env.τ, u"s"),
             au2si_round(env.toff, u"s"), au2si_round(env.tmax, u"s"))

function TruncatedGaussianEnvelope(field_params::Dict{Symbol,Any}, carrier)
    test_field_parameters(field_params, [:T]) # Period time required to round time window up
    test_field_parameters(field_params, [:τ, :σ, :α])
    test_field_parameters(field_params, [:toff, :σoff])
    test_field_parameters(field_params, [:σmax, :tmax, :Tmax])

    gaussian_common!(field_params, carrier, Tmax_rounder=NoUnits)
    @namespace!(field_params) do
        if σoff
            σoff = austrip(σoff)
            toff = σ*σoff
        else
            toff = austrip(toff)
        end
    end

    @unpack τ, σ, α, toff, tmax, Tmax = field_params
    toff < tmax || throw(ArgumentError("Hard turn-off must occur before end of pulse"))
    TruncatedGaussianEnvelope(austrip(τ), austrip(σ), α, toff, tmax, Tmax)
end

duration(env::TruncatedGaussianEnvelope) = env.τ
continuity(::TruncatedGaussianEnvelope) = Inf # This is not exactly true
span(env::TruncatedGaussianEnvelope) = -env.tmax..env.tmax

# TODO: Take truncation into account
time_integral(env::TruncatedGaussianEnvelope) = env.σ*√(2π)

# *** Spectrum

time_bandwidth_product(::AbstractGaussianEnvelope) = 2log(2)/π

@doc raw"""
    spectrum(env::AbstractGaussianEnvelope)

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
function spectrum(env::AbstractGaussianEnvelope)
    α = env.α
    N = 1/√(2α)
    β = 1/(4α)
    ω -> N*exp(-β*ω^2)
end

# ** Trapezoidal

@doc raw"""
    TrapezoidalEnvelope

This is a very simple piecewise linear function:

```math
f(t)=
\begin{cases}
r/r_{\textrm{up}},
& 0 \leq r < r_{\textrm{up}},\\
1,
& r_{\textrm{up}} \leq r < r_{\textrm{up}} + r_{\textrm{flat}}, \\
1 - \frac{r-r_{\textrm{up}}-r_{\textrm{flat}}}{r_{\textrm{down}}},
& r_{\textrm{up}} + r_{\textrm{flat}} \leq r < r_{\textrm{up}} + r_{\textrm{flat}} + r_{\textrm{down}}, \\
0,
& \textrm{else},
\end{cases}
\quad
r \defd t/T,
```
where ``T`` is the period time of the carrier.

Required parameters:
- `λ|T|f|ν|ω|ħω`,
- `ramp | (ramp_up & ramp_down)`,
- `flat`,
- `env=:trapezoidal | env=:tophat`.

Beware that this envelope can introduce artifacts at the ends of the
pulse, such that the electric field is non-vanishing, depending of
e.g. the phase of the carrier.
"""
struct TrapezoidalEnvelope{T} <: AbstractEnvelope
    ramp_up::T
    flat::T
    ramp_down::T
    period::T

    function TrapezoidalEnvelope(ramp_up::T, flat::T, ramp_down::T, period::T) where T
        ramp_up >= 0 || error("Negative up-ramp not supported")
        flat >= 0 || error("Negative flat region not supported")
        ramp_down >= 0 || error("Negative down-ramp not supported")
        ramp_up + flat + ramp_down > 0 || error("Pulse length must be non-zero")

        new{T}(ramp_up, flat, ramp_down, period)
    end
end
envelope_types[:trapezoidal] = TrapezoidalEnvelope
# Common alias
envelope_types[:tophat] = TrapezoidalEnvelope

function (env::TrapezoidalEnvelope{T})(t) where T
    t = real(t)/env.period
    if t < 0
        zero(T)
    elseif t < env.ramp_up
        t/env.ramp_up
    elseif t < env.ramp_up + env.flat
        one(T)
    elseif t < env.ramp_up + env.flat + env.ramp_down
        one(T) - (t-env.ramp_up-env.flat)/env.ramp_down
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
    TrapezoidalEnvelope(ramp_up, flat, ramp_down, austrip(T))
end

continuity(::TrapezoidalEnvelope) = 0
span(env::TrapezoidalEnvelope{T}) where T =
    zero(T)..((env.ramp_up + env.flat + env.ramp_down)*env.period)

function duration(env::TrapezoidalEnvelope)
    s = span(env)
    s.right-s.left
end

# This is formally correct, but not very useful
time_bandwidth_product(::TrapezoidalEnvelope) = Inf

# ** Cos²

@doc raw"""
    Cos²Envelope

```math
f(t) = \begin{cases}
\cos^2 \left(\frac{\pi t}{cT}\right),
& -1 \leq t/(cT) \leq 1,\\
0, & \textrm{else},
\end{cases}
```
where ``c`` is the number of cycles from zero to zero of the
``\cos^2`` envelope and ``T`` the period time.

Required parameters:
- `λ|T|f|ν|ω|ħω`,
- `cycles`,
- `env=:cos² | env=:cos2`.
"""
struct Cos²Envelope{T} <: AbstractEnvelope
    cycles::T
    period::T

    function Cos²Envelope(cycles::T, period::T) where T
        cycles >= 0 || error("Negative duration not supported")
        new{T}(cycles, period)
    end
end
envelope_types[:cos²] = Cos²Envelope
envelope_types[:cos2] = Cos²Envelope

function (env::Cos²Envelope{T})(t) where T
    t /= (env.cycles*env.period)
    if -1 ≤ 2real(t) ≤ 1
        complex_squared(cospi(t))
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
    Cos²Envelope(cycles, austrip(T))
end

duration(env::Cos²Envelope) = env.cycles*env.period
continuity(::Cos²Envelope) = 0
function span(env::Cos²Envelope)
    s = env.cycles*env.period/2
    -s..s
end

# This is formally correct, but not very useful
time_bandwidth_product(::Cos²Envelope) = Inf

# ** Exports

export duration, continuity
