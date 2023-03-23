# * Analytic Fourier transforms

function spectrum end
function convolution end

struct DiracComb{T,U}
    frequencies::Vector{Tuple{T,U}}
end

@doc raw"""
    convolution(f̂::Function, dc::DiracComb, ω)

Evaluate the convolution between the function `f̂` (assumed to be the
Fourier transform of a function `f`) with the [`DiracComb`](@ref)
`dc`. This is used to implement the Fourier transform of a function
product ``f(t)g(t)``, where ``g(t)`` is a sum of monochromatic waves:

```math
f(t)g(t)
\rightsquigarrow
\frac{1}{\sqrt{2\pi}}
(\hat{f}\star\hat{g})(\omega)
```
"""
function convolution(f::Function, dc::DiracComb{T,U}, ω) where {T,U}
    v = zero(U)

    η = inv(√(2T(π)))
    for (ω₀,c) in dc.frequencies
        v += η * c * f(ω-ω₀)
    end

    v
end

# * FFT

fftω(args...) = 2π*fftfreq(args...)
fftω(t::AbstractRange) = fftω(length(t), 1/step(t))

FFTW.fft(f::AbstractField, t::AbstractRange) =
    fft(field_amplitude(f, t))
fft_vector_potential(f::AbstractField, t::AbstractRange) =
    fft(vector_potential(f, t))

function _nfft(Y::AbstractVector, ts...)
    f = one(eltype(Y))
    for t in ts
        f *= (t[end]-t[1])/(length(t)*√(2π))
    end

    f*Y
end

@doc raw"""
    nfft(y, ts...)

Normalized FFT of `y` with respect to the axes `ts`. We use the
symmetric normalization ``(2\pi)^{-N/2}`` where ``N`` is the number of
dimensions.
"""
nfft(y::AbstractVector, ts...) = _nfft(fft(y), ts...)

nfft(f::AbstractField, t::AbstractRange) = _nfft(fft(f, t), t)

nfft_vector_potential(f::AbstractField, t::AbstractRange) =
    _nfft(fft_vector_potential(f, t), t)

# * Exports

export spectrum,
    fftω,
    fft, fft_vector_potential,
    nfft, nfft_vector_potential
