abstract type AbstractField end

Base.parent(f::AbstractField) = f

abstract type Polarization end
struct LinearPolarization <: Polarization end
struct ArbitraryPolarization <: Polarization end
Base.Broadcast.broadcastable(p::Polarization) = Ref(p)

"""
    AbstractCarrier

Any carrier, of any dimensionality.
"""
abstract type AbstractCarrier end
"""
    LinearCarrier <: AbstractCarrier

Any carrier which is one-dimensional, i.e. linearly polarized (by
convention along ``z``).
"""
abstract type LinearCarrier <: AbstractCarrier end
"""
    TransverseCarrier <: AbstractCarrier

Any carrier which is two-dimensional, i.e. polarized in the plane
perpendicular to the direction of propagation.
"""
abstract type TransverseCarrier <: AbstractCarrier end

abstract type AbstractEnvelope end

Base.Broadcast.broadcastable(c::AbstractCarrier) = Ref(c)
Base.Broadcast.broadcastable(e::AbstractEnvelope) = Ref(e)
Base.Broadcast.broadcastable(f::AbstractField) = Ref(f)

vector_potential(::LinearPolarization, f, t::AbstractVector) =
    vector_potential.(f, t)
vector_potential(::ArbitraryPolarization, f, t::AbstractVector) =
    transpose(reduce(hcat, vector_potential.(f, t)))
vector_potential(f::AbstractField, t::AbstractVector) =
    vector_potential(polarization(f), f, t)

@doc raw"""
    vector_potential(f, t)

Compute the vector potential ``\vec{A}(t)`` of the field `f`.
"""
function vector_potential end

field_amplitude(::LinearPolarization, f, t::AbstractVector) =
    field_amplitude.(f, t)
field_amplitude(::ArbitraryPolarization, f, t::AbstractVector) =
    transpose(reduce(hcat, field_amplitude.(f, t)))
field_amplitude(f::AbstractField, t::AbstractVector) =
    field_amplitude(polarization(f), f, t)

intensity(f::AbstractField, t::AbstractVector) =
    intensity.(f, t)

complex_derivative(f::Function, t::Real) = ForwardDiff.derivative(f, t)

function complex_derivative(f::Function, z::Complex)
    # https://discourse.julialang.org/t/automatic-differentiation-of-complex-valued-functions/30263/3
    ff = ((x,y),) -> begin
        fz = f(complex(x,y))
        [real(fz), imag(fz)]
    end
    J = ForwardDiff.jacobian(ff, [real(z), imag(z)])
    # Complex derivative for an analytic function
    #   f(x + im*y) = u(x,y) + im*v(x,y)
    # is given by
    #   f′(x + im*y) = uₓ + im*vₓ = v_y - im*u_y
    J[1,1] + im*J[2,1]
end

@doc raw"""
    field_amplitude(f, t)

Compute the field amplitude from the [`vector_potential`](@ref)
``A(t)`` using automatic differentiation according to

```math
F(t) = -\partial_t A(t).
```

"""
field_amplitude(f::AbstractField, t) =
    -complex_derivative(Base.Fix1(vector_potential, f), t)

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

@doc raw"""
    instantaneous_intensity(f, t)

Compute the instantaneous intensity of the field `f` as

```math
I_i(t) = \abs{F(t)}^2 = \abs{-\partial_t A(t)}^2.
```
"""
instantaneous_intensity(f::AbstractField, t::Number) = norm(field_amplitude(f, t))^2

function intensity(::LinearPolarization, f, t::Number; kwargs...)
    fun = ϕ -> -instantaneous_intensity(phase_shift(f, ϕ), t)
    res = optimize(fun, 0.0, 2π; kwargs...)
    -Optim.minimum(res)
end

function intensity(::ArbitraryPolarization, f, t::Number; kwargs...)
    cf = d -> begin
        fun = ϕ -> -abs2(field_amplitude(phase_shift(f, ϕ), t)[d])
        res = optimize(fun, 0.0, 2π; kwargs...)
        -Optim.minimum(res)
    end
    sum(cf, 1:3)
end

"""
    intensity(f, t)

Compute the intensity _envelope_ of the field `f` at time `t` by
applying [`phase_shift`](@ref) to the [`carrier`](@ref) and looking
for the maximal [`instantaneous_intensity`](@ref).
"""
intensity(f::AbstractField, t::Number) = intensity(polarization(f), f, t)

"""
    field_envelope(f, t)

Compute the field amplitude envelope as the square root of
[`field_envelope`](@ref).
"""
field_envelope(f::AbstractField, t::Number) = √(intensity(f, t))

@doc raw"""
    wavelength(f)

Return the carrier wavelength ``\lambda`` of `f`.
"""
wavelength(f::AbstractField) = wavelength(carrier(f))

"""
    period(f)

Return the carrier period time ``T`` of `f`.
"""
period(f::AbstractField) = period(carrier(f))

"""
    frequency(f)

Return the carrier frequency ``f`` of the field `f`.
"""
frequency(f::AbstractField) = frequency(carrier(f))

"""
    max_frequency(f)

Return the maximum carrier frequency ``f`` in case of a composite
field `f`.
"""
max_frequency(f::AbstractField) = max_frequency(carrier(f))

@doc raw"""
    wavenumber(f)

Return the carrier wavenumber ``\nu`` of the field `f`.
"""
wavenumber(f::AbstractField) = wavenumber(carrier(f))
fundamental(f::AbstractField) = fundamental(carrier(f))

@doc raw"""
    photon_energy(f)

Return the carrier photon energy ``\hbar\omega`` of the field `f`.
"""
photon_energy(f::AbstractField) = photon_energy(carrier(f))

@doc raw"""
    duration(f)

Return the pulse duration ``\tau`` of the field `f`.
"""
duration(f::AbstractField) = duration(envelope(f))

@doc raw"""
    continuity(f)

Return the pulse continuity, i.e. differentiability, of the field `f`.
"""
continuity(f::AbstractField) = continuity(envelope(f))

@doc raw"""
    fluence(f)

Compute the fluence of the field `f`, i.e.
```math
\frac{1}{\hbar\omega}
\int\diff{t} I(t),
```
where ``I(t)`` is the [`intensity`](@ref) envelope of the pulse.
"""
function fluence(F::AbstractField)
    ∫t = time_integral(F)
    Iau = austrip(3.5094452e16*u"W"/(u"cm"^2))

    ∫t*Iau*intensity(F)/photon_energy(F)
end

"""
    intensity(f)

Return the peak intensity of the field `f`.
"""
function intensity end

"""
    amplitude(f)

Return the peak amplitude of the field `f`.
"""
function amplitude end

"""
    carrier(f)

Return the carrier of the field `f`.
"""
function carrier end

"""
    envelope(f)

Return the envelope of the field `f`.
"""
function envelope end

"""
    phase(c)

Return the phase of the carrier `c`.
"""
function phase end

"""
    phase_shift(f, δϕ)

Return a new field whose [`carrier`](@ref) has been phase-shifted by
`δϕ`.
"""
function phase_shift end

"""
    dimensions(f)

Return the number of dimensions of the field `f`.
"""
function dimensions end

# * Linear field

"""
    LinearField

Linearly polarized field, i.e. the field amplitude is
scalar. Consisting of a peak vector potential, a vector envelope, and
a carrier, which has to be a [`LinearCarrier`](@ref).

# Examples

```jldoctest
julia> @field(A) do
       I₀ = 1e14u"W/cm^2"
       λ = 800.0u"nm"
       τ = 6.2u"fs"
       tmax = 20.0u"fs"
       end
Linearly polarized field with
  - I₀ = 2.8495e-03 au = 1.0e14 W cm⁻² =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m⁻¹
    - A₀ = 0.9372 au
  – a Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV)
  – and a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)

julia> @field(B) do
       I₀ = 0.05
       ω = 1.0
       ramp = 1.0
       flat = 3.0
       env = :trapezoidal
       end
Linearly polarized field with
  - I₀ = 5.0000e-02 au = 1.7547226e15 W cm⁻² =>
    - E₀ = 2.2361e-01 au = 114.9832 GV m⁻¹
    - A₀ = 0.2236 au
  – a Fixed carrier @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV)
  – and a /1‾3‾1\\ cycles trapezoidal envelope
```
"""
struct LinearField{Carrier<:LinearCarrier,Envelope,T} <: AbstractField
    carrier::Carrier
    env::Envelope # Vector potential envelope
    I₀::T
    E₀::T
    A₀::T
    params::Dict{Symbol, Any}
end

function LinearField(carrier, env, params)
    @unpack I₀, E₀, A₀ = params

    LinearField(carrier, env, Iaustrip(I₀), austrip(E₀), austrip(A₀), params)
end

vector_potential(f::LinearField, t::Number) = f.A₀*f.env(t)*f.carrier(t)

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
             f.E₀, au2si_round(f.E₀, u"V/m"),
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

rotation_matrix(f::LinearField{<:Any,<:Any,T}) where T = SMatrix{3,3,T}(I)

# * Transverse field

@doc raw"""
    TransverseField

Transversely polarized field, i.e the polarization vector is confined
to the plane that is perpendicular to the direction of propagation. If
the propagation vector is not rotated, this plane is the ``z-x``, with
``z`` by convention being the principal polarization axis, and ``y``
is the direction of propagation.

The formal definition is
```math
\vec{A}(t)\defd A_0 f(t) \Im\{\vec{J}\exp[\im(\omega t + \phi)]\},
```
using the complex-valued [Jones vector](https://en.wikipedia.org/wiki/Jones_calculus#Jones_vector) ``\vec{J}``,
but in practice, we use a rotation matrix ``\mat{R}`` and the [Carriers](@ref)
give the amplitudes in the canonical, unrotated coordinate sytem:
```math
\vec{A}(t) = A_0 f(t) \mat{R}\vec{C}(t).
```

# Examples

```jldoctest
julia> # Circularly polarized field

julia> @field(A) do
       I₀ = 1e14u"W/cm^2"
       λ = 800.0u"nm"
       τ = 6.2u"fs"
       tmax = 20.0u"fs"
       ξ = 1.0
       end
Transversely polarized field with
  - I₀ = 2.8495e-03 au = 1.0e14 W cm⁻² =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m⁻¹
    - A₀ = 0.9372 au
  – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV)
  – and a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)

julia> # Linearly polarized field, but explicitly in 3D

julia> @field(B) do
       I₀ = 1e14u"W/cm^2"
       λ = 800.0u"nm"
       τ = 6.2u"fs"
       tmax = 20.0u"fs"
       kind = :transverse
       end
Transversely polarized field with
  - I₀ = 2.8495e-03 au = 1.0e14 W cm⁻² =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m⁻¹
    - A₀ = 0.9372 au
  – a LinearTransverseCarrier: Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV)
  – and a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)

julia> # Linearly polarized field, rotated

julia> @field(C) do
       I₀ = 1e14u"W/cm^2"
       λ = 800.0u"nm"
       τ = 6.2u"fs"
       tmax = 20.0u"fs"
       rotation = π/3, [0,0,1]
       end
Transversely polarized field with
  - I₀ = 2.8495e-03 au = 1.0e14 W cm⁻² =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m⁻¹
    - A₀ = 0.9372 au
  – a LinearTransverseCarrier: Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV)
  – a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)
  – and a rotation of 0.33π about [0.000, 0.000, 1.000]
```
"""
struct TransverseField{Carrier<:TransverseCarrier,Envelope,Rotation,T} <: AbstractField
    carrier::Carrier
    env::Envelope # Vector potential envelope
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

vector_potential(f::TransverseField, t::Number) = f.A₀*f.env(t)*(f.R*f.carrier(t))

intensity(f::TransverseField) = f.I₀
amplitude(f::TransverseField) = f.E₀

rotation_matrix(f::TransverseField) = f.R

function show(io::IO, f::TransverseField)
    printfmt(io, """
Transversely polarized field with
  - I₀ = {1:.4e} au = {2:s} =>
    - E₀ = {3:.4e} au = {4:s}
    - A₀ = {5:.4f} au
  – a """,
             f.I₀, Iau*f.I₀,
             f.E₀, au2si_round(f.E₀, u"V/m"),
             f.A₀)
    show(io, f.carrier)
    isrotated = f.R isa AbstractMatrix
    if isrotated
        write(io, " \n  – a ")
        show(io, f.env)
        printfmt(io, " \n  – and a rotation of {1:.2f}π about [{2:.3f}, {3:.3f}, {4:.3f}]",
                 rotation_angle(f.R)/π, rotation_axis(f.R)...)
    else
        write(io, " \n  – and a ")
        show(io, f.env)
    end
end

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

# * Constant field

@doc raw"""
    ConstantField(tmax, E₀)

The field amplitude of a constant field is defined as
```math
F(t) = \begin{cases}
E_0, & 0 \leq t \leq t_{\textrm{max}}, \\
0, & \textrm{else},
\end{cases}
\implies
A(t) = \begin{cases}
-E_0t, & 0 \leq t \leq t_{\textrm{max}}, \\
0, & \textrm{else}.
\end{cases}
```

Since the vector potential is non-zero at the end of the pulse, this
is a _non-propagating_ field, i.e. it does not correspond to a freely
propagating pulse. It however corresponds to the field in an idealized
capacitor, i.e. two plates of opposite charge.

# Example

```jldoctest
julia> @field(F) do
       I₀ = 1e13u"W/cm^2"
       tmax = 3.0u"fs"
       kind = :constant
       end
Constant field of
  - 124.02412000600552 jiffies = 3.0000 fs duration, and
  - E₀ = 1.6880e-02 au = 8.680210982018475e9 V m⁻¹
```
"""
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

field_amplitude(f::ConstantField, t::Number) =
    f.E₀*(0 ≤ t && t ≤ f.tmax)

intensity(f::ConstantField, t) = field_amplitude(f, t)^2

function vector_potential(f::ConstantField{T}, t::Number) where T
    t = clamp(t, zero(T), f.tmax)
    -f.E₀*t
end

intensity(f::ConstantField) = f.E₀^2
amplitude(f::ConstantField) = f.E₀

polarization(::ConstantField) = LinearPolarization()

duration(f::ConstantField) = f.tmax
span(f::ConstantField{T}) where T = zero(T)..f.tmax

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
    field_envelope,
    params, dimensions,
    phase_shift, phase
