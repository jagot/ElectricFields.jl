# * Analytic Fourier transforms

function spectrum end
function convolution end

struct DiracComb{T,U}
    frequencies::Vector{Tuple{T,U}}
end

function Base.:(*)(R::RT, dc::DiracComb{T,U}) where {RT<:Union{SMatrix{3,3},UniformScaling},T,U<:SVector{3}}
    frequencies = Vector{Tuple{T,U}}(undef, length(dc.frequencies))
    for (i,(ω₀,c)) in enumerate(dc.frequencies)
        frequencies[i] = (ω₀, R*c)
    end
    DiracComb(frequencies)
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
"""
    fftω(t::AbstractRange)

Return the angular frequency grid corresponding to uniform sampling in
time with the tempolar grid `t`.
"""
fftω(t::AbstractRange) = fftω(length(t), 1/step(t))

"""
    fft(f::AbstractField, t::AbstractRange)

Compute the FFT of the [`field_amplitude`](@ref) of the field `f`
sampled on the uniform temporal grid `t`.
"""
FFTW.fft(f::AbstractField, t::AbstractRange) =
    fft(field_amplitude(f, t), 1)

"""
    fft_vector_potential(f::AbstractField, t::AbstractRange)

Compute the FFT of the [`vector_potential`](@ref) of the field `f`
sampled on the uniform temporal grid `t`.
"""
fft_vector_potential(f::AbstractField, t::AbstractRange) =
    fft(vector_potential(f, t), 1)

function _nfft(Y, ts...)
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
nfft(y, ts...) = _nfft(fft(y, 1:length(ts)), ts...)

"""
    nfft(f::AbstractField, t::AbstractRange)

Compute the symmetrically normalized FFT of the
[`field_amplitude`](@ref) of the field `f` sampled on the uniform
temporal grid `t`.
"""
nfft(f::AbstractField, t::AbstractRange) = _nfft(fft(f, t), t)

"""
    nfft_vector_potential(f::AbstractField, t::AbstractRange)

Compute the symmetrically normalized FFT of the
[`vector_potential`](@ref) of the field `f` sampled on the uniform
temporal grid `t`.
"""
nfft_vector_potential(f::AbstractField, t::AbstractRange) =
    _nfft(fft_vector_potential(f, t), t)

# * Exports

export spectrum,
    fftω, fftshift,
    fft, fft_vector_potential,
    nfft, nfft_vector_potential
