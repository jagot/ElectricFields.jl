# Dispersion

Dispersion is most easily calculated in the frequency domain, where it
simply amounts to multiplication of the spectral amplitudes by a
transfer function:
```math
\hat{\vec{F}}_2(\omega) =
H(\omega)
\hat{\vec{F}}_1(\omega)
```
where ``H(\omega)`` is the transfer function, ``\hat{\vec{F}}_1(\omega)`` is
the field before dispersion, and ``\hat{\vec{F}}_2(\omega)`` after. Since
``\hat{\vec{F}}(\omega)\equiv-\im\omega\hat{\vec{A}}(\omega)``, the same relation
holds for the vector potential:
```math
\hat{\vec{A}}_2(\omega) =
H(\omega)
\hat{\vec{A}}_1(\omega).
```

Dispersion is implemented in ElectricFields.jl via the transform pair
[`rfft`](@ref)/[`irfft`](@ref):

```math
\vec{A}_1\{t\}
\overset{\texttt{rfft}}{\to}
\hat{\vec{A}}_1\{\omega\}
\to
\hat{\vec{A}}_2\{\omega\} =
H\{\omega\}
\hat{\vec{A}}_1\{\omega\}
\overset{\texttt{irfft}}{\to}
\vec{A}_2\{t\}
```
where the notation ``\{t\}`` and ``\{\omega\}`` signifies uniformly
spaced vectors, since the FFT is a numeric algorithm. Since
[`field_amplitude`](@ref) and [`vector_potential`](@ref) allow the
evaluation of a field at arbitrary time points (and integrals over
arbitrary time intervals), we fit ``\vec{A}_2\{t\}`` to a
[`BSplineField`](@ref). The next complication is that the
transfer function ``H`` will in general introduce a temporal spread,
such that the compact support of the dispersed field is much wider
than that of the original one. We therefore successively increase the
time span (using [`ElectricFields.find_time_span`](@ref)) we evaluate
``\vec{A}_2\{t\}`` on until it has converged. Before we fit this
function to a [`BSplineField`](@ref), we truncate the time span again to a window
where ``|\vec{A}_2\{t\}|>\epsilon``.

```@docs
DispersedField
```

## Simple dispersive elements

Any purely dispersive medium (i.e. no loss or gain) can be written as
```math
H(\omega) = \exp[-\im\phi(\omega)]
```
where the phase can be Taylor-expanded as
```math
\phi(\omega) =
\phi_0 +
a(\omega-\omega_0) +
b(\omega-\omega_0)^2 +
...,
```
``\phi_0`` is a [`PhaseShift`](@ref), ``a=2\pi t_d`` amounts to a time
[`delay`](@ref) by ``t_d`` (by the Fourier shift theorem), and ``b``
introduced a [`Chirp`](@ref).

```julia-repl
julia> @field(F) do
           λ = 800u"nm"
           I₀ = 1.0
           τ = 3u"fs"
           σoff = 4.0
           σmax = 6.0
           env = :trunc_gauss
           ϕ = π
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
    - A₀ = 17.5580 au
  – a Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz); CEP = 1.00π
  – and a Truncated Gaussian envelope of duration 124.0241 jiffies = 3.0000 fs (intensity FWHM; turn-off from 5.0959 fs to 7.6439 fs)
  – and a bandwidth of 0.0224 Ha = 608.3170 meV ⟺ 147.0904 THz ⟺ 5933.9307 Bohr = 314.0101 nm
  – Uₚ = 77.0706 Ha = 2.0972 keV => α = 308.2823 Bohr = 16.3136 nm

julia> Fc = chirp(F, austrip(5u"fs^2"), verbosity=4)
┌ Info: Finding large enough time span to encompass dispersed field
│   f =
│    Linearly polarized field with
│      - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
│        - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
│        - A₀ = 17.5580 au
│      – a Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz); CEP = 1.00π
│      – and a Truncated Gaussian envelope of duration 124.0241 jiffies = 3.0000 fs (intensity FWHM; turn-off from 5.0959 fs to 7.6439 fs)
│      – and a bandwidth of 0.0224 Ha = 608.3170 meV ⟺ 147.0904 THz ⟺ 5933.9307 Bohr = 314.0101 nm
│      – Uₚ = 77.0706 Ha = 2.0972 keV => α = 308.2823 Bohr = 16.3136 nm
│   de = Chirp(b = 8545.5457 = 5.0000 fs², ω₀ = 0.0570 = 1.5498 eV)
│   max_iter = 7
│   ξ = 2.0
└   tol = 0.0005
----------------------------------------------------------------------------------------------------
1
t′ = -632.0183332958276:1.1049271561115868:633.1232604519392
R = 2.095947834009663
----------------------------------------------------------------------------------------------------
2
t′ = -1264.0366665916551:1.1049271561115868:1265.1415937477668
R = 0.16138039790711892
----------------------------------------------------------------------------------------------------
3
t′ = -2528.0733331833103:1.1049271561115868:2529.178260339422
R = 0.0013929399886955529
----------------------------------------------------------------------------------------------------
4
t′ = -5056.146666366621:1.1049271561115868:5057.251593522733
R = 0.00012089884366390845
┌ Info: Truncated to time interval -896.38155595379 .. 888.6605640985183
│   a = 1480
│   b = 3098
│   cutoff = 0.0014901161193847656
└   abs_cutoff = 0.014115901307391281
┌ Info: Generated B-spline
│   num_knots = 324
└   B = BSpline basis with typename(ElectricFields.LinearKnotSet)(Float64) of order k = 3 (parabolic) on -896.38155595379 .. 888.6605640985183 (324 intervals)
DispersedField:
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
    - A₀ = 17.5580 au
  – a Fixed carrier @ λ = 800.0000 nm (T = 2.6685 fs, ω = 0.0570 Ha = 1.5498 eV, f = 374.7406 THz); CEP = 1.00π
  – and a Truncated Gaussian envelope of duration 124.0241 jiffies = 3.0000 fs (intensity FWHM; turn-off from 5.0959 fs to 7.6439 fs)
  – and a bandwidth of 0.0224 Ha = 608.3170 meV ⟺ 147.0904 THz ⟺ 5933.9307 Bohr = 314.0101 nm
  – Uₚ = 77.0706 Ha = 2.0972 keV => α = 308.2823 Bohr = 16.3136 nm
  – dispersed through Chirp(b = 8545.5457 = 5.0000 fs², ω₀ = 0.0570 = 1.5498 eV)
```

![Chirped field](figures/chirped_field.svg)

In black, we see the original field, in red the chirped field, and the
bottom panel shows a [Gabor
transform](https://en.wikipedia.org/wiki/Gabor_transform) of the
chirped field, along with horizontal line at the carrier frequency
``\omega_0``, and a diagonal line at the expected instantaneous
frequency ``\omega=\omega_0 + \frac{1}{2}\frac{b}{\gamma^2 + b^2}t``,
where ``\gamma = \frac{\tau^2}{8\ln2}``.

```@docs
ElectricFields.DispersiveElement
PhaseShift
phase_shift(f::ElectricFields.AbstractField, ϕ; kwargs...)
Chirp
chirp
ElectricFields.CascadedDispersiveElement
ElectricFields.find_time_span
```

## Dispersive media

[Index ellipsoid](https://en.wikipedia.org/wiki/Index_ellipsoid)

[Interactive demonstration](https://micro.magnet.fsu.edu/primer/java/polarizedlight/ellipsoid/index.html)

```@docs
ElectricFields.Medium
IsotropicMedium
IsotropicMedium(material, d; ω₀=nothing)
Crystal
Crystal(material, d, R=I; ω₀=nothing)
```

### Sellmeier equations

```@docs
ElectricFields.Sellmeier
BK7
SiO₂
Calcite
Quartz
KTP
```

## B-splines

The present B-spline implementation is cannibalized from
[CompactBases.jl](https://github.com/JuliaApproximation/CompactBases.jl),
to avoid dragging in all the dependencies.

```@docs
BSplineField
ElectricFields.BSpline
```
