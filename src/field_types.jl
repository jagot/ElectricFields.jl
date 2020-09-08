abstract type AbstractField end
abstract type AbstractCarrier end
abstract type AbstractEnvelope end

Base.Broadcast.broadcastable(c::AbstractCarrier) = Ref(c)
Base.Broadcast.broadcastable(e::AbstractEnvelope) = Ref(e)
Base.Broadcast.broadcastable(f::AbstractField) = Ref(f)

(f::AbstractField)(t::Number) = envelope(f)(t)*carrier(f)(t)
(f::AbstractField)(t::AbstractVector{<:Number}) = f.(t)

vector_potential(f::AbstractField, t) = vector_potential(envelope(f), t)*carrier(f)(t)

@doc raw"""
    field_amplitude(f, t)

Compute the field amplitude from the [`vector_potential`](@ref)
``A(t)`` using automatic differentiation according to

```math
F(t) = -\partial_t A(t).
```

"""
field_amplitude(f::AbstractField, t) =
    -first(Zygote.gradient(Base.Fix1(vector_potential, f), t))

@doc raw"""
    field_amplitude(f, a, b)

Compute the time integral of the field amplitude according to

```math
\int_a^b\diff{t}F(t) = -[A(b) - A(a)],
```

where ``A(t)`` is the [`vector_potential`](@ref).
"""
field_amplitude(f::AbstractField, a, b) =
    -(vector_potential(f, b) - vector_potential(f, a))

wavelength(f::AbstractField) = wavelength(carrier(f))
period(f::AbstractField) = period(carrier(f))

frequency(f::AbstractField) = frequency(carrier(f)) |> base_units[:f]
max_frequency(f::AbstractField) = max_frequency(carrier(f)) |> base_units[:f]
wavenumber(f::AbstractField) = wavenumber(carrier(f))
fundamental(f::AbstractField) = fundamental(carrier(f))
photon_energy(f::AbstractField) = photon_energy(carrier(f)) |> base_units[:ħω]


intensity(f::AbstractField) = intensity(envelope(f))
amplitude(f::AbstractField) = amplitude(envelope(f))
duration(f::AbstractField) = duration(envelope(f))
continuity(f::AbstractField) = continuity(envelope(f))

# * Linear field

struct LinearField <: AbstractField
    carrier::AbstractCarrier
    env::AbstractEnvelope # Amplitude envelope
    params::Dict{Symbol, Any}
end

(f::LinearField)(t::Unitful.Time) = f.carrier(t)*f.env(t)
(f::LinearField)(fs::Unitful.Frequency=default_sampling_frequency(f)) =
    f(timeaxis(f, fs))

function show(io::IO, f::LinearField)
    write(io, "Linearly polarized field with\n  – a ")
    show(io, f.carrier)
    write(io, " \n  – and a ")
    show(io, f.env)
end

carrier(f::LinearField) = f.carrier
envelope(f::LinearField) = f.env
params(f::LinearField) = f.params
dimensions(::LinearField) = 2

function phase_shift(f::LinearField, δϕ)
    carrier = phase_shift(f.carrier, δϕ)
    p = copy(params(f))
    p[:ϕ] = phase(carrier)
    LinearField(carrier, f.env, p)
end

# * Transverse field

struct TransverseField <: AbstractField
    z::LinearField
    x::LinearField
end

duration(f::TransverseField) = max(duration.((f.z,f.x))...)
dimensions(::TransverseField) = 3

# * Keldysh parameter

@doc raw"""
    keldysh(f, Iₚ)

The [Keldysh
parameter](https://en.wikipedia.org/wiki/Tunnel_ionization) relates
the strength of a dynamic electric field to that of the binding
potential of an atom. It is given by

```math
\gamma = \sqrt{\frac{I_p}{2U_p}},
```

where ``I_p`` is the ionization potential of the atom and ``U_p`` is
the ponderomotive potential of the dynamic field.
"""
keldysh(f::AbstractField, Iₚ::Unitful.Energy) =
    √(Iₚ/2params(f)[:Uₚ]) |> NoUnits

# * Exports

export carrier,
    wavelength, period,
    frequency, max_frequency, wavenumber, fundamental, photon_energy,
    envelope,
    intensity, amplitude,
    duration,
    field_amplitude, vector_potential,
    params, dimensions,
    phase_shift, phase,
    keldysh
