function spectrum end
function convolution end

struct DiracComb{T,U}
    frequencies::Vector{Tuple{T,U}}
end

function convolution(f::Function, dc::DiracComb{T,U}, ω) where {T,U}
    v = zero(U)

    for (ω₀,c) in dc.frequencies
        v += c * f(ω-ω₀)
    end

    v
end

export spectrum
