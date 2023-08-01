```@meta
CurrentModule = ElectricFields
DocTestSetup = quote
    using ElectricFields, Unitful, UnitfulAtomic
end
```

# ElectricFields.jl

The idea of this package is to provide an interface between the
reality and calculations. In the calculations, it is useful to
represent fields in terms of cycles of a _fundamental_ frequency,
which yields a timebase. E.g. one might use laser pulses of 800 nm, in
an experiment, which has a period time of about 2.66 fs. It is,
however, easier to calculate in normalized time, and relate all other
quantities of interest (such as ionization potential, &c) to this time
scale.

The package provides a simple DSL that requires just a handful of
parameters, that can be given in _any_ unit system (thanks to
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl)). Different
fields can be combined in any way that is physically reasonable, to
recreate complicated experimental situations. Everything can then be
converted to normalized time, for use in calculations.

## Scope

At the moment, only electric fields that are the solution of the
[Helmholtz
equation](https://en.wikipedia.org/wiki/Helmholtz_equation),
i.e. those that are separable in the time and space coordinates, are
supported. Further more, the spatial dependence is neglected entirely
(long wavelength approximation), and only the temporal behaviour is
modelled. This appropriate when studying systems which are smaller
than the fundamental wavelength of the field, such as atoms in a laser
field. It may be problematic for larger molecules, which can be
sensitive to the spatial variation of the field.

The temporal structure of the field is modelled as
```math
\vec{A}(t) = A_0 f(t) \Im\{\vec{J} \exp[\im(\omega t + \phi)]\},
```
where ``A_0`` is the peak _vector potential_, ``f(t)`` is the
_envelope_ (which does not have to be slowly varying in time, i.e. can
be ultrashort), ``\vec{J}`` is the [_Jones vector_](https://en.wikipedia.org/wiki/Jones_calculus#Jones_vector), and
``\exp[\im(\omega t + \phi)]`` is the _carrier_. At present only
linear polarization, or polarization transverse to the direction of
propagation is supported, i.e. no longitudinal polarization is
implemented.

The electric field amplitude is computed as
```math
\vec{F}(t) \defd -\partial_t \vec{A}(t)
```
via automatic differentiation
([ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)).

Presently, only the temporal picture is well supported in
ElectricFields.jl, since the main usage (for now) is to provide input
to time-propagation of the time-dependent Schrödinger equation. The
equivalent spectral picture, related to the temporal one via the
Fourier transform, will receive better support in a future release.

## Examples

### Linear polarization

The simplest pulse is created thus:

```jldoctest index_examples
julia> using ElectricFields, Unitful

julia> @field(IR) do
       I₀ = 1e14u"W/cm^2"
       λ = 800.0u"nm"
       τ = 6.2u"fs"
       σmax = 6.0
       end
Linearly polarized field with
  - I₀ = 2.8495e-03 au = 1.0e14 W cm^-2 =>
    - E₀ = 5.3380e-02 au = 27.4492 GV m^-1
    - A₀ = 0.9372 au
  – a Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz)
  – and a Gaussian envelope of duration 6.2000 fs (intensity FWHM; ±6.08σ)
  – and a bandwidth of 0.0108 Ha = 294.3469 meV ⟺ 71.1728 THz ⟺ 2871.2568 Bohr = 151.9404 nm
  – Uₚ = 0.2196 Ha = 5.9759 eV => α = 16.4562 Bohr = 870.8242 pm

julia> vector_potential(IR, 4.0)
0.21160647961322301

julia> field_amplitude(IR, 4.0), field_envelope(IR, 4.0)
(-0.05194703530281701, 0.05336224984361336)

julia> instantaneous_intensity(IR, 4.0), intensity(IR, 4.0)
(0.0026984944767521166, 0.002847529708372214)

julia> span(IR)
-661.9198939608042 .. 661.9198939608042

julia> timeaxis(IR)
-661.9198939608042:1.1041199232040104:661.9198939608042
```

![Simple example](figures/index_example.svg)

If we give no units, [Hartree atomic
units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are assumed:
```jldoctest index_examples
julia> @field(XUV) do
       I₀ = 0.04
       ω = 1.0
       τ = 150.0
       σmax = 4
       end
Linearly polarized field with
  - I₀ = 4.0000e-02 au = 1.40377808e15 W cm^-2 =>
    - E₀ = 2.0000e-01 au = 102.8441 GV m^-1
    - A₀ = 0.2000 au
  – a Fixed carrier @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV, f = 6.5797 PHz)
  – and a Gaussian envelope of duration 3.6283 fs (intensity FWHM; ±4.04σ)
  – and a bandwidth of 0.0185 Ha = 502.9732 meV ⟺ 121.6184 THz ⟺ 15.9151 Bohr = 842.1896 pm
  – Uₚ = 0.0100 Ha = 272.1138 meV => α = 0.2000 Bohr = 10.5835 pm
```

### Arbitrary polarization, composite fields

Other [Envelopes](@ref) and polarizations, as well as some simple arithmetic is possible:
```jldoctest index_examples
julia> @field(A) do
           I₀ = 1.0
           ω = 1.0
           cycles = 6.0
           env = :cos²
           ξ = 1.0
       end
Transversely polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 1.0000 au
  – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV, f = 6.5797 PHz)
  – and a 6.00 cycles cos² envelope
  – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
  – Uₚ = 0.2500 Ha = 6.8028 eV => α = 1.0000 Bohr = 52.9177 pm

julia> @field(B) do
           I₀ = 1.0
           ω = 2.0
           cycles = 6.0
           env = :cos²
           ξ = -1.0
       end
Transversely polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.5000 au
  – a Elliptical carrier with ξ = -1.00 (LCP) @ λ = 22.7817 nm (T = 75.9915 as, ω = 2.0000 Ha = 54.4228 eV, f = 13.1594 PHz)
  – and a 6.00 cycles cos² envelope
  – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
  – Uₚ = 0.0625 Ha = 1.7007 eV => α = 0.2500 Bohr = 13.2294 pm

julia> F = A + delay(B, 3/2π)
┌ Transversely polarized field with
│   - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
│     - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
│     - A₀ = 1.0000 au
│   – a Elliptical carrier with ξ = 1.00 (RCP) @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV, f = 6.5797 PHz)
│   – and a 6.00 cycles cos² envelope
│   – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
│   – Uₚ = 0.2500 Ha = 6.8028 eV => α = 1.0000 Bohr = 52.9177 pm
⊕
│ Transversely polarized field with
│   - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
│     - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
│     - A₀ = 0.5000 au
│   – a Elliptical carrier with ξ = -1.00 (LCP) @ λ = 22.7817 nm (T = 75.9915 as, ω = 2.0000 Ha = 54.4228 eV, f = 13.1594 PHz)
│   – and a 6.00 cycles cos² envelope
│   – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
│   – Uₚ = 0.0625 Ha = 1.7007 eV => α = 0.2500 Bohr = 13.2294 pm
└   – delayed by 0.4775 jiffies = 11.5493 as


julia> field_amplitude(F, 4.0)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.8793235934912674
 -0.0
  0.06802883592577502
```

![Simple polarized example](figures/index_polarized_example.svg)

### Spectra

We can compute the spectrum of a field using either the Fast Fourier
Transform [`fft`](@ref), [`nfft`](@ref), or for some fields
analytically [`field_amplitude_spectrum`](@ref):

```jldoctest
julia> @field(F) do
           I₀ = 1.0
           T = 1.0
           τ = 2.0
           σmax = 6.0
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
    - A₀ = 0.1592 au
  – a Fixed carrier @ λ = 7.2516 nm (T = 24.1888 as, ω = 6.2832 Ha = 170.9742 eV, f = 41.3414 PHz)
  – and a Gaussian envelope of duration 48.3777 as (intensity FWHM; ±7.06σ)
  – and a bandwidth of 1.3863 Ha = 37.7230 eV ⟺ 9.1214 PHz ⟺ 30.2350 Bohr = 1.6000 nm
  – Uₚ = 0.0063 Ha = 172.3181 meV => α = 0.0253 Bohr = 1.3404 pm

julia> t = timeaxis(F)
-6.0:0.010008340283569641:6.0

julia> Fv = field_amplitude(F, t);

julia> Av = vector_potential(F, t);

julia> ω = fftshift(fftω(t));

julia> # We need to undo the phase, since the FFT does not care that
       # pulse is centred around zero.
       F̂v = exp.(im*ω*t[1]) .* fftshift(nfft(F, t), 1);

julia> Âv = exp.(im*ω*t[1]) .* fftshift(nfft_vector_potential(F, t), 1);

julia> F̂v_exact = field_amplitude_spectrum(F, ω);

julia> Âv_exact = vector_potential_spectrum(F, ω);
```

![Simple spectrum example](figures/index_spectrum_example.svg)

## Reference

```@index
```

```@meta
DocTestSetup = nothing
```
