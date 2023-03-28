Base.parent(f::AbstractField) = f

abstract type Polarization end
struct LinearPolarization <: Polarization end
struct ArbitraryPolarization <: Polarization end
Base.Broadcast.broadcastable(p::Polarization) = Ref(p)

field_types = Dict{Symbol,Any}()

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

# Complex derivative for an analytic function
#   f(x + im*y) = u(x,y) + im*v(x,y)
# is given by
#   f′(x + im*y) = uₓ + im*vₓ = v_y - im*u_y
assemble_derivative(::Number, J) = J[1,1] + im*J[2,1]
assemble_derivative(::SVector{3}, J) = SVector{3}(J[1:3,1] + im*J[4:6,1])

function complex_derivative(f::Function, z::Complex)
    # https://discourse.julialang.org/t/automatic-differentiation-of-complex-valued-functions/30263/3
    ff = ((x,y),) -> begin
        fz = f(complex(x,y))
        vcat(real(fz), imag(fz))
    end
    J = ForwardDiff.jacobian(ff, [real(z), imag(z)])
    assemble_derivative(f(z), J) # This amounts to an extra function
                                 # call, but what can you do
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

"""
    time_bandwidth_product(f)

Return the time–bandwidth product of the field `f`.
"""
time_bandwidth_product(f::AbstractField) = time_bandwidth_product(envelope(f))

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

Return the number of dimensions of the field `f`. See also
[`polarization`](@ref).

# Examples

```jldoctest
julia> @field(F) do
           I₀ = 2.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 2.0000e+00 au = 7.0188904e16 W cm^-2 =>
    - E₀ = 1.4142e+00 au = 727.2178 GV m^-1
    - A₀ = 0.4502 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm

julia> dimensions(F)
1

julia> @field(F) do
           I₀ = 2.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
           ξ = 1.0
       end
Transversely polarized field with
  - I₀ = 2.0000e+00 au = 7.0188904e16 W cm^-2 =>
    - E₀ = 1.4142e+00 au = 727.2178 GV m^-1
    - A₀ = 0.4502 au
  – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm

julia> dimensions(F)
3
```
"""
function dimensions end

"""
    polarization(f)

Return the polarization kind of the field `f`. See also
[`dimensions`](@ref).

# Examples

```jldoctest
julia> @field(F) do
           I₀ = 2.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
       end
Linearly polarized field with
  - I₀ = 2.0000e+00 au = 7.0188904e16 W cm^-2 =>
    - E₀ = 1.4142e+00 au = 727.2178 GV m^-1
    - A₀ = 0.4502 au
  – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm

julia> polarization(F)
LinearPolarization()

julia> @field(F) do
           I₀ = 2.0
           T = 2.0
           σ = 3.0
           Tmax = 3.0
           ξ = 1.0
       end
Transversely polarized field with
  - I₀ = 2.0000e+00 au = 7.0188904e16 W cm^-2 =>
    - E₀ = 1.4142e+00 au = 727.2178 GV m^-1
    - A₀ = 0.4502 au
  – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
  – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±2.00σ)
  – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
  – Uₚ = 0.0507 Ha = 1.3785 eV => α = 0.1433 Bohr = 7.5826 pm

julia> polarization(F)
ArbitraryPolarization()
```
"""
function polarization end

function show_bandwidth(io::IO, f::AbstractField)
    τ = duration(f)
    tbp = time_bandwidth_product(f)
    λ = wavelength(f)
    ω = photon_energy(f)
    Δf = tbp/τ
    Δω = Δf*2π
    Δλ = (2π*austrip(1u"c")/ω^2)*Δω
    printfmt(io, "and a bandwidth of {1:.4f} Ha = {2:s} ⟺ {3:s} ⟺ {4:.4f} Bohr = {5:s}",
             Δω, au2si_round(Δω, u"eV"),
             au2si_round(Δf, u"Hz"),
             Δλ, au2si_round(Δλ, u"m"))
end

# ** Spectra

"""
    vector_potential_spectrum(f::AbstractField, ω)

Compute the analytic Fourier transform of the
[`vector_potential`](@ref) of the field `f` at the angular frequency
`ω`.
"""
function vector_potential_spectrum end

vector_potential_spectrum(::LinearPolarization, f, ω::AbstractVector) =
    vector_potential_spectrum.(f, ω)
vector_potential_spectrum(::ArbitraryPolarization, f, ω::AbstractVector) =
    transpose(reduce(hcat, vector_potential_spectrum.(f, ω)))
vector_potential_spectrum(f::AbstractField, ω::AbstractVector) =
    vector_potential_spectrum(polarization(f), f, ω)

@doc raw"""
    field_amplitude_spectrum(f::AbstractField, ω)

Compute the analytic Fourier transform of the
[`field_amplitude`](@ref) of the field `f` at the angular frequency
`ω`, using the Fourier identity

```math
\vec{F}(t) = -\partial_t \vec{A}(t)
\iff
\hat{\vec{F}}(\omega) = -\im\omega \hat{\vec{A}}(\omega)
```
"""
field_amplitude_spectrum(f::AbstractField, ω::Number) =
    -im*ω*vector_potential_spectrum(f, ω)

field_amplitude_spectrum(::LinearPolarization, f, ω::AbstractVector) =
    field_amplitude_spectrum.(f, ω)
field_amplitude_spectrum(::ArbitraryPolarization, f, ω::AbstractVector) =
    transpose(reduce(hcat, field_amplitude_spectrum.(f, ω)))
field_amplitude_spectrum(f::AbstractField, ω::AbstractVector) =
    field_amplitude_spectrum(polarization(f), f, ω)

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
  - I₀ = 2.8495e-03 au = 1.0e14 W cm^-2 =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m^-1
    - A₀ = 0.9372 au
  – a Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz)
  – and a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)
  – and a bandwidth of 0.0108 Ha = 294.3469 meV ⟺ 71.1728 THz ⟺ 2871.2568 Bohr = 151.9404 nm
  – Uₚ = 0.2196 Ha = 5.9759 eV => α = 16.4562 Bohr = 870.8242 pm

julia> @field(B) do
       I₀ = 0.05
       ω = 1.0
       ramp = 1.0
       flat = 3.0
       env = :trapezoidal
       end
Linearly polarized field with
  - I₀ = 5.0000e-02 au = 1.7547226e15 W cm^-2 =>
    - E₀ = 2.2361e-01 au = 114.9832 GV m^-1
    - A₀ = 0.2236 au
  – a Fixed carrier @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV, f = 6.5797 PHz)
  – and a /1‾3‾1\\ cycles trapezoidal envelope
  – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
  – Uₚ = 0.0125 Ha = 340.1423 meV => α = 0.2236 Bohr = 11.8328 pm
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
vector_potential(f::LinearField) = f.A₀

vector_potential_spectrum(f::LinearField, ω::Number) =
    f.A₀*convolution(spectrum(f.env), spectrum(f.carrier), ω)

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
    write(io, "\n  – and a ")
    show(io, f.env)
    write(io, "\n  – ")
    show_bandwidth(io, f)
    write(io, "\n  – ")
    show_strong_field_properties(io, f)
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
  - I₀ = 2.8495e-03 au = 1.0e14 W cm^-2 =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m^-1
    - A₀ = 0.9372 au
  – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz)
  – and a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)
  – and a bandwidth of 0.0108 Ha = 294.3469 meV ⟺ 71.1728 THz ⟺ 2871.2568 Bohr = 151.9404 nm
  – Uₚ = 0.2196 Ha = 5.9759 eV => α = 16.4562 Bohr = 870.8242 pm

julia> # Linearly polarized field, but explicitly in 3D

julia> @field(B) do
       I₀ = 1e14u"W/cm^2"
       λ = 800.0u"nm"
       τ = 6.2u"fs"
       tmax = 20.0u"fs"
       kind = :transverse
       end
Transversely polarized field with
  - I₀ = 2.8495e-03 au = 1.0e14 W cm^-2 =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m^-1
    - A₀ = 0.9372 au
  – a LinearTransverseCarrier: Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz)
  – and a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)
  – and a bandwidth of 0.0108 Ha = 294.3469 meV ⟺ 71.1728 THz ⟺ 2871.2568 Bohr = 151.9404 nm
  – Uₚ = 0.2196 Ha = 5.9759 eV => α = 16.4562 Bohr = 870.8242 pm

julia> # Linearly polarized field, rotated

julia> @field(C) do
       I₀ = 1e14u"W/cm^2"
       λ = 800.0u"nm"
       τ = 6.2u"fs"
       tmax = 20.0u"fs"
       rotation = π/3, [0,0,1]
       end
Transversely polarized field with
  - I₀ = 2.8495e-03 au = 1.0e14 W cm^-2 =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m^-1
    - A₀ = 0.9372 au
  – a LinearTransverseCarrier: Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz)
  – a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±8.11σ)
  – and a rotation of 0.33π about [0.000, 0.000, 1.000]
  – and a bandwidth of 0.0108 Ha = 294.3469 meV ⟺ 71.1728 THz ⟺ 2871.2568 Bohr = 151.9404 nm
  – Uₚ = 0.2196 Ha = 5.9759 eV => α = 16.4562 Bohr = 870.8242 pm
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
vector_potential(f::TransverseField) = f.A₀

vector_potential_spectrum(f::TransverseField, ω::Number) =
    f.A₀*convolution(spectrum(f.env), f.R*spectrum(f.carrier), ω)

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
        write(io, "\n  – a ")
        show(io, f.env)
        printfmt(io, "\n  – and a rotation of {1:.2f}π about [{2:.3f}, {3:.3f}, {4:.3f}]",
                 rotation_angle(f.R)/π, rotation_axis(f.R)...)
    else
        write(io, "\n  – and a ")
        show(io, f.env)
    end
    write(io, "\n  – ")
    show_bandwidth(io, f)
    write(io, "\n  – ")
    show_strong_field_properties(io, f)
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
A(t_{\textrm{max}}), & t_{\textrm{max}} < t, \\
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
  - 124.0241 jiffies = 3.0000 fs duration, and
  - E₀ = 1.6880e-02 au = 8.6802 GV m^-1
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
  - {1:.4f} jiffies = {2:s} duration, and
  - E₀ = {3:.4e} au = {4:s}""",
             f.tmax, au2si_round(f.tmax, u"s"),
             f.E₀, au2si_round(f.E₀, u"V/m"))
end

intensity(f::ConstantField, t::Number) = field_amplitude(f, t)^2

function vector_potential(f::ConstantField{T}, t::Number) where T
    z = zero(T)
    if t < z
        z
    elseif t < f.tmax
        -f.E₀*t
    else
        -f.E₀*f.tmax
    end
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

field_types[:constant] = ConstantField

function field_amplitude_spectrum(f::ConstantField, ω::Number)
    a = inv(f.tmax)
    exp(-im*a*ω/2)*(1/√(2π*a^2))*sinc(ω/(2π*a))
end

function vector_potential_spectrum(f::ConstantField, ω::Number)
    # # This is not entirely correct
    # a = inv(f.tmax)
    # b = 1/(2π*a)
    # -im*exp(-im*ω/2)*(1/√(2π*a^2))*(-im*sinc(b*ω)/2 + cosc(b*ω)*b)
    -field_amplitude_spectrum(f,ω)/(im*ω)
end

# * Ramps

@doc raw"""
    Ramp(tmax, E₀, f, name)

The field amplitude of a ramp is defined as

```math
F(t) = \begin{cases}
E_0f'(\tau), & 0 \leq t_{\textrm{max}}, \\
0, \textrm{else},
\end{cases}
\implies
A(t) = \begin{cases}
-E_0t_{\textrm{max}}f(\tau), & 0 \leq t \leq t_{\textrm{max}},\\
A(t_{\textrm{max}}), & t_{\textrm{max}} < t,\\
0, & \textrm{else},
\end{cases}
\quad
\tau = \frac{t}{t_{\textrm{max}}}.
```

To define a new ramp, one thus only needs to define one function on
the unit interval, ``f(\tau)`` whose derivative ``f'(\tau)`` rises
from ``0`` to ``1``. Similar to [`ConstantField`](@ref), `Ramp` is a
_non-propagating_ field, but is realizable in e.g. a capacitor.

Three kinds of ramps are predefined, `:linear_ramp`,
`:parabolic_ramp`, `:sin²_ramp` (with the alias `:sin2_ramp`).

# Examples

```jldoctest
julia> @field(F) do
           E₀ = 1.0
           tmax = 4.0u"fs"
           kind = :sin²_ramp
       end
sin² up-ramp of
  - 165.3655 jiffies = 4.0000 fs duration, and
  - E₀ = 1.0000e+00 au = 514.2207 GV m^-1

julia> @field(F) do
           E₀ = 1.0
           tmax = 4.0u"fs"
           kind = :linear_ramp
           ramp = :down
       end
Linear down-ramp of
  - 165.3655 jiffies = 4.0000 fs duration, and
  - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
```
"""
struct Ramp{T,AFunc} <: AbstractField
    tmax::T
    E₀::T
    f::AFunc
    name::String
    params::Dict{Symbol, Any}
end

function Ramp(f, name, field_params)
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

    r = get(field_params, :ramp, :up)
    a,suffix = if r == :up
        f, "up"
    elseif r == :down
        A = f(1)
        τ -> A - f(1-τ), "down"
    else
        throw(ArgumentError("Unknown :ramp kind $(r)"))
    end

    Ramp(austrip(tmax), austrip(E₀), a, name*" "*suffix*"-", field_params)
end

function Base.show(io::IO, f::Ramp)
    printfmt(io, """
$(f.name)ramp of
  - {1:.4f} jiffies = {2:s} duration, and
  - E₀ = {3:.4e} au = {4:s}""",
             f.tmax, au2si_round(f.tmax, u"s"),
             f.E₀, au2si_round(f.E₀, u"V/m"))
end

ElectricFields.field_types[:linear_ramp] =
    p -> Ramp(t -> t^2/2, "Linear", p)

ElectricFields.field_types[:parabolic_ramp] =
    p -> Ramp(t -> t^2 - t^3/3, "Parabolic", p)

ElectricFields.field_types[:sin2_ramp] =
    ElectricFields.field_types[:sin²_ramp] =
    p -> Ramp(t -> t/2 - sinpi(t)/2π, "sin²", p)

intensity(f::Ramp, t::Number) = field_amplitude(f, t)^2

function vector_potential(f::Ramp{T}, t::Number) where T
    z = zero(T)
    if t < z
        z
    elseif t < f.tmax
        τ = clamp(t, zero(T), f.tmax)/f.tmax
        -f.E₀*f.tmax*f.f(τ)
    else
        -f.E₀*f.tmax*f.f(one(T))
    end
end

polarization(::Ramp) = LinearPolarization()

duration(f::Ramp) = f.tmax
span(f::Ramp{T}) where T = zero(T)..f.tmax

# These methods are just for convenience (to be able to establish a
# time base for calculations), but they are not physical as such.
period(f::Ramp{T}) where T = one(T)
max_frequency(f::Ramp{T}) where T = one(T)
photon_energy(f::Ramp{T}) where T = 2T(π)

dimensions(::Ramp) = 1

# * Exports

export LinearPolarization, ArbitraryPolarization, polarization,
    carrier,
    wavelength, period,
    frequency, max_frequency, wavenumber, fundamental, photon_energy,
    envelope,
    intensity, amplitude, fluence,
    duration, time_bandwidth_product,
    field_amplitude, vector_potential, instantaneous_intensity,
    field_envelope,
    params, dimensions,
    phase_shift, phase,
    field_amplitude_spectrum, vector_potential_spectrum
