# Tabulated fields

Depending on which quantity that is tabulated, we choose different
wrapper types:

| Quantity         | Type                                        |
| ---------------- | ------------------------------------------- |
| ``F_i(q)``       | [`TabulatedFieldAmplitude`](@ref)           |
| ``\|F_i(q)\|``   | [`TabulatedFieldMagnitude`](@ref)           |
| ``\|F_i(q)\|^2`` | [`TabulatedFieldPowerSpectrum`](@ref)       |
| ``A_i(q)``       | [`TabulatedVectorPotentialAmplitude`](@ref) |
| ``\|A_i(q)\|^2`` | [`TabulatedVectorPotentialMagnitude`](@ref) |

The superscript ``i`` refers to the Cartesian components
``x,y,z``.

1. The first step is to fit the tabulated quantity to a
   [`ElectricFields.BSpline`](@ref).
2. Depending on the unit of `q` and the physical quantity provided, we then
   follow different tracks to ``A(t)`` which is what the final field
   expects, e.g.:
   ```math
   \hat{\vec{F}}\{\omega\} \to
   \hat{\vec{A}}\{\omega\} =
   -\frac{1}{\im\omega}\hat{\vec{F}}\{\omega\}
   \overset{\texttt{irfft}}{\to} \vec{A}\{t\}
   ```
   If we start out with a quantity given as a magnitude to a power
   ``n`` and a phase, we first find the complex amplitude as
   ```math
   |V_i(q)|^n,\exp[{\im\phi(q)}] \to
   V_i(q) = [|V_i(q)|^n]^{1/n}\exp[{\im\phi(q)}].
   ```
   If the parameter is not an angular frequency ``\omega``, we then
   make the appropriate variable change, via `rfft` if the signal is
   given in the time domain.
   ```math
   \begin{aligned}
   V_i\{t\} &\overset{\texttt{rfft}}{\to} \hat{V}_i\{\omega\}, &
   \hat{V}_i\{f\} &\to \hat{V}_i\{\omega\} = \hat{V}_i\left\{2\pi f\right\}, &
   \hat{V}_i\{\lambda\} &\to \hat{V}_i\{\omega\} = \hat{V}_i\left\{\frac{2\pi c}{\lambda}\right\}.
   \end{aligned}
   ```
3. ``\vec{A}\{t\}`` is subsequently fit to a new B-spline and wrapped
   in a [`BSplineField`](@ref).

```@docs
tabulated_field
ElectricFields.TabulatedFieldQuantity
TabulatedFieldAmplitude
TabulatedFieldMagnitude
TabulatedFieldPowerSpectrum
TabulatedVectorPotentialAmplitude
TabulatedVectorPotentialMagnitude
```
