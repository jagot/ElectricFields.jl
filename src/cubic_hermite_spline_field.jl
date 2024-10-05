"""
    CubicHermiteSplineField(t, A, Aₜ)

Represents an electric field using cubic Hermite spline interpolation
of known values, where `t` are the sample times, `A` the array (vector
for 1d fields, 3-column matrix for 3d fields) of sample values for the
vector potential, and `Aₜ` a similar array for minus the electric
field, since the cubic Hermite spline interpolation requires the
derivative at each sample node.
"""
struct CubicHermiteSplineField{Tt<:AbstractVector,At<:AbstractVector} <: AbstractField
    t::Tt
    A::At
    Aₜ::At
    function CubicHermiteSplineField(t::Tt, A::At, Aₜ::At) where {Tt<:AbstractVector,At<:AbstractVector}
        n = length(t)
        nA = size(A, 1)
        nAₜ = size(Aₜ, 1)
        n == nA == nAₜ ||
            throw(DimensionMismatch("Number of nodes $(n) must match the number of function values $(nA) and derivatives $(nAₜ)"))
        new{Tt,At}(t, A, Aₜ)
    end
end

function CubicHermiteSplineField(t, A::AbstractMatrix, Aₜ::AbstractMatrix)
    size(A, 2) == size(Aₜ, 2) == 3 || throw(ArgumentError("Three components expected for matrix input"))
    n = length(t)
    nA = size(A, 1)
    nAₜ = size(Aₜ, 1)
    n == nA == nAₜ ||
        throw(DimensionMismatch("Number of nodes $(n) must match the number of function values $(nA) and derivatives $(nAₜ)"))

    CubicHermiteSplineField(t, [SVector(A[i,1],A[i,2],A[i,2]) for i = 1:n],
                            [SVector(Aₜ[i,1],Aₜ[i,2],Aₜ[i,2]) for i = 1:n])
end

"""
    CubicHermiteSplineField(t, A, ::Nothing)

Convenience constructor for [`CubicHermiteSplineField`](@ref), when
the vector potential is known. Its derivative (minus the electric
field) will be approximated; if `t` is uniformly spaced, a FFT-based
derivative will computed (it is up to the user to apodize the samples,
if needed).
"""
CubicHermiteSplineField(t::AbstractVector, A::AbstractVecOrMat, ::Nothing) =
    CubicHermiteSplineField(t, A, approximate_derivative(t, A))

"""
    CubicHermiteSplineField(t, ::Nothing, F)

Convenience constructor for [`CubicHermiteSplineField`](@ref), when
the electric field is known. Its integral (minus the vector potential)
will be approximated; if `t` is uniformly spaced, a FFT-based integral
will computed (it is up to the user to apodize the samples, if
needed).
"""
CubicHermiteSplineField(t::AbstractVector, ::Nothing, F::AbstractVecOrMat) =
    CubicHermiteSplineField(t, -approximate_integral(t, F), -F)

Base.show(io::IO, f::CubicHermiteSplineField) =
    printfmt(io, "{1:d}-sample, {2:d}-component Cubic Hermite spline field, t ∈ {3:s}",
             length(f.t), dimensions(f), span(f))

function Base.show(io::IO, ::MIME"text/plain", f::CubicHermiteSplineField)
    show(io, f)
    println(io)
    print(io, "  ")
    show_strong_field_properties(io, f)
end

span(f::CubicHermiteSplineField) = first(f.t)..last(f.t)

min_step(t::AbstractRange) = step(t)
min_step(t::AbstractVector) = minimum(abs, diff(t))
max_frequency(f::CubicHermiteSplineField) = inv(min_step(f.t))

dimensions(::CubicHermiteSplineField{<:Any,<:AbstractVector}) = 1
dimensions(::CubicHermiteSplineField{<:Any,<:AbstractMatrix}) = 3

polarization(::CubicHermiteSplineField{<:Any,<:AbstractVector}) = LinearPolarization()
polarization(::CubicHermiteSplineField{<:Any,<:AbstractMatrix}) = ArbitraryPolarization()

vector_potential(f::CubicHermiteSplineField, t::Number) =
    cubic_hermite_interpolation(f.t, f.A, f.Aₜ, t)

export CubicHermiteSplineField

# * Strong field properties

function _maximize_interpolation(t, v, f, f′; verbosity=0)
    # This is only a rough approximation, and will most likely only
    # work reasonably well if there is one isolated pulse in the
    # trace.

    N = length(t)

    # We begin by finding the maximum point.
    i = argmax(i -> abs(v[i]), 1:N)
    tᵢ = t[i]
    fᵢ = abs(v[i])
    verbosity > 0 && @info "Initial maximum amplitude" i tᵢ fᵢ

    # Since the electric field is minus the time derivative of the
    # vector potential, we know that when the vector potential is
    # extremized, the electric field has a root. The converse is not
    # necessarily true (only for purely monochromatic fields is that
    # guaranteed), but we assume it is.
    tₘₐₓ = find_zero(f′, tᵢ)
    fₘₐₓ = abs(f(tₘₐₓ))

    verbosity > 0 && @info "Optimized maximum amplitude" tₘₐₓ fₘₐₓ fₘₐₓ ≥ fᵢ

    if fₘₐₓ < fᵢ
        @warn "Was not able to find the proper maximum"
        return fᵢ
    end

    fₘₐₓ
end

function amplitude(f::CubicHermiteSplineField; kwargs...)
    # Since the electric field is minus the time derivative of the
    # vector potential, we know that when the vector potential is
    # extremized, the electric field has a root. The converse is not
    # necessarily true (only for purely monochromatic fields is that
    # guaranteed), but we assume it is.
    _maximize_interpolation(f.t, f.Aₜ,
                            Base.Fix1(field_amplitude, f),
                            (Base.Fix1(vector_potential, f),
                             t -> -field_amplitude(f, t));
                            kwargs...)
end

function ponderomotive_potential(f::CubicHermiteSplineField; kwargs...)
    # Since the electric field is minus the time derivative of the
    # vector potential, we know that when the vector potential is
    # extremized, the electric field has a root. We therefore look for
    # a root of the electric field, starting from the time we found
    # above.
    Aₘₐₓ = _maximize_interpolation(f.t, f.A,
                                   Base.Fix1(vector_potential, f),
                                   t -> -field_amplitude(f, t);
                                   kwargs...)

    Aₘₐₓ^2/4
end

function free_oscillation_amplitude(f::CubicHermiteSplineField; kwargs...)
    ∫A = approximate_integral(f.t, f.A)
    ∫Af = t -> cubic_hermite_interpolation(f.t, ∫A, f.A, t)

    _maximize_interpolation(f.t, ∫A,
                            ∫Af,
                            (Base.Fix1(vector_potential, f),
                             t -> -field_amplitude(f, t));
                            kwargs...)
end

# * Cubic Hermite spline interpolation

function cubic_hermite_interpolation(xp::AbstractVector, f::AbstractVector, fₓ::AbstractVector, x::Number)
    i = find_interval(xp, x)
    i == length(xp) && return f[i]

    δx = (xp[i+1] - xp[i])
    t = (x - xp[i])/δx

    # See https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Representations
    a = 1 - t
    b = a^2
    c = t^2

    h₀₀ = (1+2t)*b
    h₁₀ = t*b
    h₀₁ = c*(3 - 2t)
    h₁₁ = -c*a

    h₀₀*f[i] + h₁₀*δx*fₓ[i] + h₀₁ * f[i+1] + h₁₁*δx*fₓ[i+1]
end

cubic_hermite_interpolation(xp::AbstractVector, f::AbstractVector, fₓ::AbstractVector, x::AbstractVector) =
    cubic_hermite_interpolation.(Ref(xp), Ref(f), Ref(fₓ), x)

# * FFT derivatives/integrals

# The derivative/integral of a real function is real, but using the
# FFT to compute these operation necessarily introduces an imaginary
# part, which ideally should be "small".
function ensure_real(y::AbstractArray{Complex{T}}; tol=√(eps(T))) where T
    reY = real(y)
    imY = imag(y)
    norm(imY)/norm(reY) < tol || @warn "Expected a negligible imaginary part, got" norm(reY) norm(imY) norm(imY)/norm(reY)
    reY
end

#=
Here we implement FFT-based differentition/integration, as described
by

- Johnson, S. G. (2011). Notes on FFT-based
  differentiation. https://math.mit.edu/~stevenj/fft-deriv.pdf
=#

function fft_derivative(y::AbstractVecOrMat{<:Complex}, fs=1)
    # This implements Algorithm 1 by Johnson (2011).
    N = size(y, 1)
    ω = 2π*fftfreq(N, fs)
    Y = fft(y, 1)
    Y′ = im*ω .* Y
    if iseven(N)
        Y′[N÷2,:] .= false
    end
    ifft(Y′, 1)
end

fft_derivative(y::AbstractVecOrMat{<:Real}, args...) =
    ensure_real(fft_derivative(complex(y), args...))

function fft_integral(y::AbstractVecOrMat{<:Complex}, fs=1)
    # This implements the inverse of Algorithm 1 by Johnson (2011).
    N = size(y, 1)
    ω = 2π*fftfreq(N, fs)
    Y = fft(y, 1)
    Y′ = Y ./ (im*ω)
    # We have to force the DC component to zero, otherwise we're
    # trying to integrate a constant over an infinite interval (due to
    # the division in the previous step, Y′[1,:] is ±∞ or NaN).
    Y′[1,:] .= false
    if iseven(N)
        Y′[N÷2,:] .= false
    end
    ifft(Y′, 1)
end

fft_integral(y::AbstractVecOrMat{<:Real}, args...) =
    ensure_real(fft_integral(complex(y), args...))

approximate_derivative(t::AbstractRange, A) = fft_derivative(A, 1/step(t))
approximate_integral(t::AbstractRange, A) = fft_integral(A, 1/step(t))
