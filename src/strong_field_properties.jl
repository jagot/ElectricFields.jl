# * Ponderomotive potential
@doc raw"""
    ponderomotive_potential(f)

Return the [ponderomotive
potential](https://en.wikipedia.org/wiki/Ponderomotive_energy)
``U_p``, which is the cycle-average quiver energy of a free electron
in an electromagnetic field `f`. It is given by

```math
U_p =
\frac{e^2E_0^2}{4m\omega^2}=\frac{2e^2}{c\varepsilon_0m}\times\frac{I}{4\omega^2},
```

or, in atomic units,

```math
U_p = \frac{I}{4\omega^2}.
```
"""
ponderomotive_potential(f::Union{LinearField,TransverseField,WrappedField}) =
    params(f)[:Uₚ]

ponderomotive_potential(f::SumField) =
    ponderomotive_potential(f.a) + ponderomotive_potential(f.b)

# * Keldysh parameter

@doc raw"""
    keldysh(f, Iₚ)

The [Keldysh
parameter](https://en.wikipedia.org/wiki/Tunnel_ionization) relates
the strength of a dynamic electric field to that of the binding
potential of an atom. It is given by

```math
\gamma = \sqrt{\frac{I_p}{2U_p}},
```

where ``I_p`` is the ionization potential of the atom and ``U_p`` is
the ponderomotive potential of the dynamic field.
"""
keldysh(f::AbstractField, Iₚ::Unitful.Energy) =
    √(Iₚ/2ponderomotive_potential(f)) |> NoUnits

# * Free oscillation amplitude

@doc raw"""
    free_oscillation_amplitude(F)

Compute the free oscillation amplitude of an electric field `F`,
i.e. the mean excursion length during one cycle of the field, defined
as

```math
\alpha \defd \frac{F}{\omega^2}
```

where `F` is the peak amplitude, i.e. this is defined for one cycle of
a monochrome field.
"""
free_oscillation_amplitude(F::AbstractField) =
    amplitude(F)/photon_energy(F)^2

free_oscillation_amplitude(F::Union{WrappedField,NegatedField,DelayedField}) =
    free_oscillation_amplitude(parent(F))

# This is just a heuristic, i.e. if `F.a` and `F.b` are of the same
# frequency, but exactly out of phase, the free oscillation amplitude
# is actually zero.
free_oscillation_amplitude(F::SumField) =
    free_oscillation_amplitude(F.a) +
    free_oscillation_amplitude(F.b)

# * Exports

export ponderomotive_potential, keldysh, free_oscillation_amplitude
