# * Carriers

carrier_types = Dict{Symbol,Any}()

max_frequency(carrier::AbstractCarrier) = frequency(carrier)

# ** Fixed carrier
#    The carrier is fixed in the sense that the instantaneous frequency
#    is constant throughout the pulse.

struct FixedCarrier{Λ,Tt,Ω,Φ} <: LinearCarrier
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
    printfmt(io, "Fixed carrier @ λ = {1:s} (T = {2:s}, ω = {3:.4f} Ha = {4:s})",
             si_round(u"m"(carrier.λ)),
             si_round(u"s"(carrier.T)),
             carrier.ω, au2si_round(carrier.ω, u"eV"))
    !iszero(carrier.ϕ) &&
        printfmt(io, "; CEP = {1:0.2f}π", carrier.ϕ/π)
end

# ** Transverse carriers

# *** Linear

struct LinearTransverseCarrier{Carrier<:LinearCarrier} <: TransverseCarrier
    carrier::Carrier
end

LinearTransverseCarrier(field_params::Dict{Symbol,Any}) =
    LinearTransverseCarrier(FixedCarrier(field_params))

function show(io::IO, carrier::LinearTransverseCarrier)
    write(io, "LinearTransverseCarrier: ")
    show(io, carrier.carrier)
end

(carrier::LinearTransverseCarrier)(t) = SVector(0, 0, carrier.carrier(t))

carrier_types[:linear] = LinearTransverseCarrier

for fun in [:wavelength, :period, :frequency, :wavenumber, :fundamental, :photon_energy, :phase]
    @eval $fun(carrier::LinearTransverseCarrier) = $fun(carrier.carrier)
end

phase_shift(c::LinearTransverseCarrier, δϕ) =
    LinearTransverseCarrier(phase_shift(c.carrier, δϕ))

# *** Elliptical

struct EllipticalCarrier{Λ,Tt,Ω,Φ,Ξ} <: TransverseCarrier
    λ::Λ
    T::Tt
    ω::Ω
    ϕ::Φ # Carrier–envelope phase, in radians
    ξ::Ξ # Amount of ellipticity, minor/major axis ratio
end

function EllipticalCarrier(field_params::Dict{Symbol,Any})
    @unpack λ, T, ω = field_params
    ϕ = get(field_params, :ϕ, 0)
    ξ = get(field_params, :ξ, 0)
    EllipticalCarrier(λ, T, austrip(ω), ϕ, ξ)
end

function show(io::IO, carrier::EllipticalCarrier)
    ξ = carrier.ξ
    pol = if iszero(ξ)
        ""
    elseif ξ == 1
        " (RCP)"
    elseif ξ == -1
        " (LCP)"
    elseif ξ > 0
        " (right)"
    else
        " (left)"
    end
    printfmt(io, "Elliptical carrier with ξ = {1:.2f}{2:s} @ λ = {3:s} (T = {4:s}, ω = {5:.4f} Ha = {6:s})",
             ξ, pol,
             si_round(u"m"(carrier.λ)),
             si_round(u"s"(carrier.T)),
             carrier.ω, au2si_round(carrier.ω, u"eV"))
    !iszero(carrier.ϕ) &&
        printfmt(io, "; CEP = {1:0.2f}π", carrier.ϕ/π)
end

function (carrier::EllipticalCarrier)(t)
    ξ = carrier.ξ
    ϕ = carrier.ω*t + carrier.ϕ
    s,c = sincos(ϕ)
    SVector(ξ*c, 0, s)/√(1 + ξ^2)
end

carrier_types[:elliptical] = EllipticalCarrier

wavelength(carrier::EllipticalCarrier) = carrier.λ
period(carrier::EllipticalCarrier) = carrier.T

frequency(carrier::EllipticalCarrier) = 1/carrier.T
wavenumber(carrier::EllipticalCarrier) = 1/carrier.λ
fundamental(carrier::EllipticalCarrier) = carrier.ω
photon_energy(carrier::EllipticalCarrier) = carrier.ω
phase(carrier::EllipticalCarrier) = carrier.ϕ

phase_shift(c::EllipticalCarrier, δϕ) =
    EllipticalCarrier(c.λ, c.T, c.ω, c.ϕ+δϕ, c.ξ)
