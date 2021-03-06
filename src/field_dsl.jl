# * Field creation
# ** Parameter calculation

@doc raw"""
    calc_params!(field_params)

This function performs the calculation of different quantities from
the information provided.

The [ponderomotive
potential](https://en.wikipedia.org/wiki/Ponderomotive_energy) ``U_p``
is the cycle-average quiver energy of a free electron in an
electromagnetic field. It is given by

```math
U_p =
\frac{e^2E_0^2}{4m\omega^2}=\frac{2e^2}{c\varepsilon_0m}\times\frac{I}{4\omega^2},
```

or, in atomic units,

```math
U_p = \frac{I}{4\omega^2}.
```
"""
function calc_params!(field_params::Dict{Symbol,Any})
    test_field_parameters(field_params, [:λ, :T, :f, :ν, :ω, :ħω])
    test_field_parameters(field_params, [:I₀, :E₀, :Uₚ])

    for k in keys(field_params)
        field_params[k] = get_unitful_quantity(field_params, k)
    end

    @namespace!(field_params) do
        if λ || T
            if λ
                T = λ/u"c"
            elseif T
                λ = T*u"c"
            end
            ν = 1/λ
            f = 1/T
            ω = 2π*u"rad"*f
            ħω = u"ħ"*ω
        else # ∝ Frequency specified
            if f || ν
                if f
                    ν = f/u"c"
                elseif ν
                    f = ν*u"c"
                end
                ω = 2π*u"rad"*f
            elseif ω || ħω
                if ω
                    ħω = u"ħ"*ω
                elseif ħω
                    ω = ħω/u"ħ"
                end
                f = ω/(2π*u"rad")
                ν = f/u"c"
            end
            T = 1/f
            λ = 1/ν
        end

        if I₀ || Uₚ
            if Uₚ
                I₀ = Uₚ / (2*u"q"^2/(u"c"*u"ε0"*u"me")) * 4ω^2
            end
            E₀ = √(2I₀/(u"ε0"*u"c"))
        elseif E₀
            I₀ = u"ε0"*u"c"/2*E₀^2
        end
        if !Uₚ
            Uₚ = 2*u"q"^2/(u"c"*u"ε0"*u"me") * I₀/4ω^2
        end
    end

    field_params
end

# ** Frontend macro

function make_field(field_params::Dict{Symbol,Any})
    calc_params!(field_params)

    # Maybe these two blocks can be implicitly deduced from the passed
    # parameters? E.g. if a chirp parameter is given, the carrier type
    # should autmatically be resolved as ChirpedCarrier. Similarly, if
    # ramp and flat are given, a trapezoidal pulse is requested.

    carrier_sym = get(field_params, :carrier,
                      :q ∉ keys(field_params) ? :fixed : :harmonic)
    carrier_sym ∉ keys(carrier_types) &&
        error("Unknown carrier type $(carrier_sym), valid choices are $(keys(carrier_types))")
    :q ∈ keys(field_params) && carrier_sym != :harmonic &&
        error("Invalid carrier type, $(carrier_sym), for field with harmonic components")
    carrier = carrier_types[carrier_sym](field_params)

    env_sym = get(field_params, :env, :gauss)
    env_sym ∉ keys(envelope_types) &&
        error("Unknown envelope type $(env_sym), valid choices are $(keys(envelope_types))")
    env = envelope_types[env_sym](field_params)

    :ξ in keys(field_params) &&
        error("Elliptical (transverse) fields not yet supported!")

    LinearField(carrier, env, field_params)
end

macro field(spec, var)
    spec.head == :-> ||
        error("Expected a block with parameters for definition of the field")
    block = spec.args[2]
    block.head == :block ||
        error("Expected a block with parameters for definition of the field")

    field_params = gensym()
    local tree = walk(spec, field_params, true)
    quote
        $field_params = Dict{Symbol,Any}()
        $(tree)()
        $(esc(var)) = make_field($field_params)
    end
end

export @field
