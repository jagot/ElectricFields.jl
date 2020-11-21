abstract type AbstractField end

abstract type Polarization end
struct LinearPolarization <: Polarization end
struct ArbitraryPolarization <: Polarization end

abstract type AbstractCarrier end
abstract type LinearCarrier <: AbstractCarrier end
abstract type TransverseCarrier <: AbstractCarrier end

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

instantaneous_intensity(f::AbstractField, t) = norm(field_amplitude(f, t))^2

field_envelope(f::AbstractField, t) = √(intensity(f, t))

wavelength(f::AbstractField) = wavelength(carrier(f))
period(f::AbstractField) = period(carrier(f))

frequency(f::AbstractField) = frequency(carrier(f))
max_frequency(f::AbstractField) = max_frequency(carrier(f))
wavenumber(f::AbstractField) = wavenumber(carrier(f))
fundamental(f::AbstractField) = fundamental(carrier(f))
photon_energy(f::AbstractField) = photon_energy(carrier(f))

duration(f::AbstractField) = duration(envelope(f))
continuity(f::AbstractField) = continuity(envelope(f))

function fluence(F::AbstractField)
    ∫t = time_integral(F)
    Iau = austrip(3.5094452e16*u"W"/(u"cm"^2))

    ∫t*Iau*intensity(F)/photon_energy(F)
end

# * Linear field

struct LinearField{Carrier<:LinearCarrier,Envelope,T} <: AbstractField
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

intensity(f::LinearField) = f.I₀
amplitude(f::LinearField) = f.E₀

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

polarization(::LinearField) = LinearPolarization()

carrier(f::LinearField) = f.carrier
envelope(f::LinearField) = f.env
params(f::LinearField) = f.params
dimensions(::LinearField) = 1

function phase_shift(f::LinearField, δϕ)
    carrier = phase_shift(f.carrier, δϕ)
    p = copy(params(f))
    p[:ϕ] = phase(carrier)
    LinearField(carrier, f.env, f.I₀, f.E₀, f.A₀, p)
end

time_integral(f::LinearField) = time_integral(envelope(f))

make_temp_field(carrier::LinearCarrier, env, params) =
    LinearField(carrier, env, 1, 1, 1, params)

function intensity(f::LinearField, t; kwargs...)
    fun = ϕ -> -instantaneous_intensity(phase_shift(f, ϕ), t)
    res = optimize(fun, 0.0, 2π; kwargs...)
    -Optim.minimum(res)
end

# * Transverse field

struct TransverseField{Carrier<:TransverseCarrier,Envelope,Rotation,T} <: AbstractField
    carrier::Carrier
    env::Envelope # Amplitude envelope
    I₀::T
    E₀::T
    A₀::T
    R::Rotation
    params::Dict{Symbol, Any}
end

function TransverseField(carrier, env, params)
    @unpack I₀, E₀, A₀, R = params

    TransverseField(carrier, env, Iaustrip(I₀), austrip(E₀), austrip(A₀), R, params)
end

vector_potential(f::TransverseField, t) = f.A₀*f.env(t)*(f.R*f.carrier(t))

intensity(f::TransverseField) = f.I₀
amplitude(f::TransverseField) = f.E₀

polarization(::TransverseField) = ArbitraryPolarization()

carrier(f::TransverseField) = f.carrier
envelope(f::TransverseField) = f.env
params(f::TransverseField) = f.params
dimensions(::TransverseField) = 3

function phase_shift(f::TransverseField, δϕ)
    carrier = phase_shift(f.carrier, δϕ)
    p = copy(params(f))
    p[:ϕ] = phase(carrier)
    TransverseField(carrier, f.env, f.I₀, f.E₀, f.A₀, f.R, p)
end

time_integral(f::TransverseField) = time_integral(envelope(f))

make_temp_field(carrier::TransverseCarrier, env, params,) =
    TransverseField(carrier, env, 1, 1, 1, I, params)

function intensity(f::TransverseField, t; kwargs...)
    cf = d -> begin
        fun = ϕ -> -abs2(field_amplitude(phase_shift(f, ϕ), t)[d])
        res = optimize(fun, 0.0, 2π; kwargs...)
        -Optim.minimum(res)
    end
    sum(cf, 1:3)
end

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
    -f.E₀*t
end

intensity(f::ConstantField) = f.E₀^2
amplitude(f::ConstantField) = f.E₀

polarization(::ConstantField) = LinearPolarization()

duration(f::ConstantField) = f.tmax
span(f::ConstantField{T}) where T = (zero(T), f.tmax)

# These methods are just for convenience (to be able to establish a
# time base for calculations), but they are not physical as such.
period(f::ConstantField{T}) where T = one(T)
max_frequency(f::ConstantField{T}) where T = one(T)
photon_energy(f::ConstantField{T}) where T = 2T(π)

dimensions(::ConstantField) = 1

# * Exports

export LinearPolarization, ArbitraryPolarization, polarization,
    carrier,
    wavelength, period,
    frequency, max_frequency, wavenumber, fundamental, photon_energy,
    envelope,
    intensity, amplitude, fluence,
    duration,
    field_amplitude, vector_potential, instantaneous_intensity,
    params, dimensions,
    phase_shift, phase
