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

export spectrum
