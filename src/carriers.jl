# * Carriers

carrier_types = Dict{Symbol,Any}()

max_frequency(carrier::AbstractCarrier) = frequency(carrier)

# ** Fixed carrier
#    The carrier is fixed in the sense that the instantaneous frequency
#    is constant throughout the pulse.

struct FixedCarrier <: AbstractCarrier
    λ::Unitful.Length
    T::Unitful.Time
    ω::Unitful.Frequency
    ϕ::Number # Carrier–envelope phase, in radians
end

(carrier::FixedCarrier)(t::Unitful.Time) = sin(carrier.ω*t + carrier.ϕ)

carrier_types[:fixed] = FixedCarrier

wavelength(carrier::FixedCarrier) = carrier.λ
period(carrier::FixedCarrier) = carrier.T

frequency(carrier::FixedCarrier) = 1/carrier.T
wavenumber(carrier::FixedCarrier) = 1/carrier.λ
fundamental(carrier::FixedCarrier) = carrier.ω
photon_energy(carrier::FixedCarrier) = carrier.ω * u"ħ" |> base_units[:ħω]
phase(carrier::FixedCarrier) = carrier.ϕ

phase_shift(c::FixedCarrier, δϕ) =
    FixedCarrier(c.λ, c.T, c.ω, c.ϕ+δϕ)

function FixedCarrier(field_params::Dict{Symbol,Any})
    @unpack λ, T, ω = field_params
    ϕ = get(field_params, :ϕ, 0)
    FixedCarrier(λ, T, ω, ϕ)
end

function show(io::IO, carrier::FixedCarrier)
    write(io, @sprintf("Fixed carrier @ λ = %0.2f %s (T = %0.2f %s)",
                       usplit(carrier.λ)..., usplit(carrier.T)...))
    if carrier.ϕ != 0
        write(io, @sprintf("; CEP = %0.2fπ", carrier.ϕ/π))
    end
end

# ** Harmonic carrier

struct HarmonicCarrier <: AbstractCarrier
    λ::Unitful.Length
    T::Unitful.Time
    ω::Unitful.Frequency
    ϕ::Number # Carrier–envelope phase, in radians
    q::AbstractVector{Int}
end

harmonics(carrier::HarmonicCarrier) = carrier.q
export harmonics

(carrier::HarmonicCarrier)(t::Unitful.Time) = sum(sin(q*(carrier.ω*t + carrier.ϕ))
                                                  for q ∈ harmonics(carrier))

carrier_types[:harmonic] = HarmonicCarrier

wavelength(carrier::HarmonicCarrier) = carrier.λ
period(carrier::HarmonicCarrier) = carrier.T

frequency(carrier::HarmonicCarrier) = 1/carrier.T
max_frequency(carrier::HarmonicCarrier) = maximum(harmonics(carrier))*frequency(carrier)
wavenumber(carrier::HarmonicCarrier) = 1/carrier.λ
fundamental(carrier::HarmonicCarrier) = carrier.ω
photon_energy(carrier::HarmonicCarrier) = carrier.ω * u"ħ" |> base_units[:ħω]

function HarmonicCarrier(field_params::Dict{Symbol,Any})
    @unpack λ, T, ω, q = field_params
    ϕ = get(field_params, :ϕ, 0)
    HarmonicCarrier(λ, T, ω, ϕ, q)
end

function show(io::IO, carrier::HarmonicCarrier)
    write(io, @sprintf("Harmonic carrier @ λ = %0.2f %s (T = %0.2f %s); q ∈ %s",
                       usplit(carrier.λ)..., usplit(carrier.T)...,
                       string(harmonics(carrier))))
    if carrier.ϕ != 0
        write(io, @sprintf("; CEP = %0.2fπ", carrier.ϕ/π))
    end
end

# ** Dispersed carriers [0/2]
# *** TODO Chirped carrier
# *** TODO Sellmeier equations
# ** Constant carrier
#    The carrier is constant in the sense that wavelength is infinite
#    and there is no oscillation, but the type still fulfils the carrier
#    interface, such that it can be used to establish a time-base, etc.

struct ConstantCarrier <: AbstractCarrier
    λ::Unitful.Length
    T::Unitful.Time
    ω::Unitful.Frequency
end

(carrier::ConstantCarrier)(t::Unitful.Time) = 1

carrier_types[:constant] = ConstantCarrier

wavelength(carrier::ConstantCarrier) = carrier.λ
period(carrier::ConstantCarrier) = carrier.T

frequency(carrier::ConstantCarrier) = 1/carrier.T
wavenumber(carrier::ConstantCarrier) = 1/carrier.λ
fundamental(carrier::ConstantCarrier) = carrier.ω
photon_energy(carrier::ConstantCarrier) = carrier.ω * u"ħ" |> base_units[:ħω]

function ConstantCarrier(field_params::Dict{Symbol,Any})
    @unpack λ, T, ω = field_params
    ConstantCarrier(λ, T, ω)
end

function show(io::IO, carrier::ConstantCarrier)
    write(io, @sprintf("Constant carrier @ λ = %0.2f %s (T = %0.2f %s)",
                       usplit(carrier.λ)..., usplit(carrier.T)...))
end
