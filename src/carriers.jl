# * Carriers

carrier_types = Dict{Symbol,Any}()

max_frequency(carrier::AbstractCarrier) = frequency(carrier)

# ** Fixed carrier
#    The carrier is fixed in the sense that the instantaneous frequency
#    is constant throughout the pulse.

struct FixedCarrier{Λ,Tt,Ω,Φ} <: AbstractCarrier
    λ::Λ
    T::Tt
    ω::Ω
    ϕ::Φ # Carrier–envelope phase, in radians
end

(carrier::FixedCarrier)(t) = sin(carrier.ω*t + carrier.ϕ)

carrier_types[:fixed] = FixedCarrier

wavelength(carrier::FixedCarrier) = carrier.λ
period(carrier::FixedCarrier) = carrier.T

frequency(carrier::FixedCarrier) = 1/carrier.T
wavenumber(carrier::FixedCarrier) = 1/carrier.λ
fundamental(carrier::FixedCarrier) = carrier.ω
photon_energy(carrier::FixedCarrier) = carrier.ω
phase(carrier::FixedCarrier) = carrier.ϕ

phase_shift(c::FixedCarrier, δϕ) =
    FixedCarrier(c.λ, c.T, c.ω, c.ϕ+δϕ)

function FixedCarrier(field_params::Dict{Symbol,Any})
    @unpack λ, T, ω = field_params
    ϕ = get(field_params, :ϕ, 0)
    FixedCarrier(λ, T, austrip(ω), ϕ)
end

function show(io::IO, carrier::FixedCarrier)
    printfmt(io, "Fixed carrier @ λ = {1:s} (T = {2:s})",
             si_round(carrier.λ),
             si_round(carrier.T))
    !iszero(carrier.ϕ) &&
        printfmt(io, "; CEP = {1:0.2f}π", carrier.ϕ/π)
end
