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
