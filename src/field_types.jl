struct LinearField <: AbstractField
    λ::Number
    T::Number
    ω::Number
    τ::Number
end

wavelength(f::LinearField) = f.λ
period(f::LinearField) = f.T

frequency(f::LinearField) = 1/f.T
wavenumber(f::LinearField) = 1/f.λ
fundamental(f::LinearField) = f.ω
energy(f::LinearField) = f.ω * u"hbar"

duration(f::LinearField) = f.τ

intensity(f::LinearField) = 0
amplitude(f::LinearField) = 0

struct TransverseField <: AbstractField
    z::LinearField
    x::LinearField
end

duration(f::TransverseField) = max(duration.((f.z,f.x))...)
