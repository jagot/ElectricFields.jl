```@meta
CurrentModule = ElectricFields
```

# ElectricFields

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

```@index
```

```@autodocs
Modules = [ElectricFields]
```
