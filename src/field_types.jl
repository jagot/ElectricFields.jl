abstract type AbstractField end
abstract type AbstractCarrier end
abstract type AbstractEnvelope end

Base.Broadcast.broadcastable(c::AbstractCarrier) = Ref(c)
Base.Broadcast.broadcastable(e::AbstractEnvelope) = Ref(e)
Base.Broadcast.broadcastable(f::AbstractField) = Ref(f)

@doc raw"""
    field_amplitude(f, t)

Compute the field amplitude from the [`vector_potential`](@ref)
``A(t)`` using automatic differentiation according to

```math
F(t) = -\partial_t A(t).
```

"""
field_amplitude(f::AbstractField, t) =
    -ForwardDiff.derivative(Base.Fix1(vector_potential, f), t)

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

instantaneous_intensity(f::AbstractField, t) = field_amplitude(f, t)^2
function intensity(f::AbstractField, t; kwargs...)
    fun = ϕ -> -instantaneous_intensity(phase_shift(f, ϕ), t)
    res = optimize(fun, 0.0, 2π; kwargs...)
    -Optim.minimum(res)
end

field_envelope(f::AbstractField, t) = √(intensity(f, t))

wavelength(f::AbstractField) = wavelength(carrier(f))
period(f::AbstractField) = period(carrier(f))

frequency(f::AbstractField) = frequency(carrier(f))
max_frequency(f::AbstractField) = max_frequency(carrier(f))
wavenumber(f::AbstractField) = wavenumber(carrier(f))
fundamental(f::AbstractField) = fundamental(carrier(f))
photon_energy(f::AbstractField) = photon_energy(carrier(f))

intensity(f::AbstractField) = intensity(envelope(f))
amplitude(f::AbstractField) = amplitude(envelope(f))
duration(f::AbstractField) = duration(envelope(f))
continuity(f::AbstractField) = continuity(envelope(f))

# * Linear field

struct LinearField{Carrier,Envelope,T} <: AbstractField
    carrier::Carrier
    env::Envelope # Amplitude envelope
    I₀::T
    E₀::T
    A₀::T
    params::Dict{Symbol, Any}
end

function LinearField(carrier, env, params)
    @unpack I₀, E₀, A₀ = params

    LinearField(carrier, env, Iaustrip(I₀), austrip(E₀), austrip(A₀), params)
end

vector_potential(f::LinearField, t) = f.A₀*f.env(t)*f.carrier(t)

function show(io::IO, f::LinearField)
    printfmt(io, """
Linearly polarized field with
  - I₀ = {1:.4e} au = {2:s} =>
    - E₀ = {3:.4e} au = {4:s}
    - A₀ = {5:.4f} au
  – a """,
             f.I₀, Iau*f.I₀,
             f.E₀, au2si(f.E₀, u"V/m"),
             f.A₀)
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
    LinearField(carrier, f.env, f.I₀, f.E₀, f.A₀, p)
end

# * Transverse field

struct TransverseField <: AbstractField
    z::LinearField
    x::LinearField
end

duration(f::TransverseField) = max(duration.((f.z,f.x))...)
dimensions(::TransverseField) = 3

# * Constant field

struct ConstantField{T} <: AbstractField
    tmax::T
    E₀::T
    params::Dict{Symbol, Any}
end

function ConstantField(field_params)
    test_field_parameters(field_params, [:I₀, :E₀])
    test_field_parameters(field_params, [:tmax])

    for k in keys(field_params)
        field_params[k] = get_unitful_quantity(field_params, k)
    end

    @namespace!(field_params) do
        if I₀
            E₀ = √(2I₀/(u"ε0"*u"c"))
        end
    end

    @unpack tmax, E₀ = field_params

    ConstantField(austrip(tmax), austrip(E₀), field_params)
end

function show(io::IO, f::ConstantField)
    printfmt(io, """
Constant field of
  - {1:s} jiffies = {2:s} duration, and
  - E₀ = {3:.4e} au = {4:s}""",
             f.tmax, au2si_round(f.tmax, u"s"),
             f.E₀, au2si(f.E₀, u"V/m"))
end

field_amplitude(f::ConstantField, t) =
    f.E₀*(0 ≤ t && t ≤ f.tmax)

intensity(f::ConstantField, t) = field_amplitude(f, t)^2

function vector_potential(f::ConstantField{T}, t) where T
    t = clamp(t, zero(T), f.tmax)
    -f.E₀/2*t^2
end

duration(f::ConstantField) = f.tmax
span(f::ConstantField{T}) where T = (zero(T), f.tmax)

period(f::ConstantField{T}) where T = one(T)
max_frequency(f::ConstantField{T}) where T = one(T)

dimensions(::ConstantField) = 1

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
    field_amplitude, vector_potential, instantaneous_intensity,
    params, dimensions,
    phase_shift, phase,
    keldysh
