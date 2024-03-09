# Field types

## Basic field types
```@docs
ElectricFields.LinearField
ElectricFields.TransverseField
ElectricFields.LinearTransverseField
ElectricFields.ConstantField
ElectricFields.Ramp
```

## Arithmetic
```@docs
ElectricFields.SumField
ElectricFields.WrappedField
ElectricFields.NegatedField
ElectricFields.DelayedField
ElectricFields.PaddedField
ElectricFields.WindowedField
ElectricFields.phase_shift
```

## Apodizing Windows

This is a generalization of [`ElectricFields.WindowedField`](@ref),
where the field is multiplied by an [_apodizing
window_](https://en.wikipedia.org/wiki/Window_function), that tapers
of towards the boundaries, and typically (but not always) is zero at
the boundaries. The window functions are non-zero on the interval
``[-1/2,1/2]``, so if other intervals are desired, they will need to
be rescaled. This happens automatically.

Since the apodizing window is multiplied with the vector potential,
the field amplitude gains an extra term due to the chain rule:
```math
\vec{A}_w(t) =
W(t)\vec{A}(t)
\implies
\vec{F}_w(t) =
-\partial_t
\vec{A}_w(t) =
W(t)\vec{F}(t) -
W'(t)\vec{A}(t).
```

### Usage

```julia-repl
julia> @field(F) do
           ω = 1.0
           I₀ = 1.0
           ramp = 0.0
           flat = 3.0
           env = :tophat
       end
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
    - A₀ = 1.0000 au
  – a Fixed carrier @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV, f = 6.5797 PHz)
  – and a /0‾3‾0\ cycles trapezoidal envelope
  – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
  – Uₚ = 0.2500 Ha = 6.8028 eV => α = 1.0000 Bohr = 52.9177 pm

julia> Fw = ApodizedField(F, 1.0, 14.0, ElectricFields.Kaiser(3))
Kaiser(α = 3) window from 1.0000 jiffies = 24.1888 as to 14.0000 jiffies = 338.6438 as of
Linearly polarized field with
  - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
    - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
    - A₀ = 1.0000 au
  – a Fixed carrier @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV, f = 6.5797 PHz)
  – and a /0‾3‾0\ cycles trapezoidal envelope
  – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
  – Uₚ = 0.2500 Ha = 6.8028 eV => α = 1.0000 Bohr = 52.9177 pm
```

![Apodized field](figures/apodized_field.svg)

```@docs
ApodizedField
span(::ApodizedField)
```

### Windows

See the [Wikipedia page on apodizing
windows](https://en.wikipedia.org/wiki/Window_function) for the formal
definitions of the windows and their frequency responses.

![Apodizing windows](figures/windows.svg)

```@docs
ElectricFields.AbstractWindow
ElectricFields.Rect
ElectricFields.Hann
ElectricFields.Hamming
ElectricFields.Blackman
ElectricFields.BlackmanExact
ElectricFields.Nuttall
ElectricFields.BlackmanNuttall
ElectricFields.BlackmanHarris
ElectricFields.Kaiser
```
