# Interface

## Creation

```@docs
ElectricFields.@field
```

## Time-dependent quantities

```@docs
ElectricFields.field_amplitude
ElectricFields.vector_potential
ElectricFields.instantaneous_intensity
ElectricFields.intensity(::LinearPolarization, f, t::Number; kwargs...)
ElectricFields.field_envelope
ElectricFields.span
ElectricFields.timeaxis
ElectricFields.steps
```

## Spectra

```@docs
ElectricFields.spectrum
ElectricFields.field_amplitude_spectrum
ElectricFields.vector_potential_spectrum
ElectricFields.AbstractFFTs.fft
ElectricFields.nfft
ElectricFields.fft_vector_potential
ElectricFields.nfft_vector_potential
ElectricFields.fftω
```
