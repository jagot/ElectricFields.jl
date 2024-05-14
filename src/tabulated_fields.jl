similar_fill(A::AbstractArray, f) = fill!(similar(A), f)

"""
    TabulatedFieldQuantity{T,N} <: AbstractArray{T,N}

Base type for all possible quantities one may tabulate for a field,
wraps around an ordinary `AbstractArray`.
"""
abstract type TabulatedFieldQuantity{T,N} <: AbstractArray{T,N} end

Base.parent(a::TabulatedFieldQuantity) = a.a
for f in [:size, :getindex, :setindex!]
    @eval Base.$f(a::TabulatedFieldQuantity, args...) = $f(parent(a), args...)
end

@doc raw"""
    TabulatedFieldAmplitude(a)

Indicates that the values contained in the array `a` are of the kind
``F(q)``.
"""
struct TabulatedFieldAmplitude{T,N,A<:AbstractArray{T,N}} <: TabulatedFieldQuantity{T,N}
    a::A
end

function to_complex_amplitude(V::TabulatedFieldAmplitude, ϕ)
    all(iszero, ϕ) ||
    @warn "Non-zero phase for tabulated amplitudes meaningless, ignoring."
    V
end

@doc raw"""
    TabulatedFieldMagnitude(a)

Indicates that the values contained in the array `a` are of the kind
``|F(q)|``.
"""
struct TabulatedFieldMagnitude{T,N,A<:AbstractArray{T,N}} <: TabulatedFieldQuantity{T,N}
    a::A
end

to_complex_amplitude(V::TabulatedFieldMagnitude, ϕ) =
    TabulatedFieldAmplitude(V.a .* exp.(im*ϕ))

@doc raw"""
    TabulatedFieldPowerSpectrum(a)

Indicates that the values contained in the array `a` are of the kind
``|F(q)|^2``.
"""
struct TabulatedFieldPowerSpectrum{T,N,A<:AbstractArray{T,N}} <: TabulatedFieldQuantity{T,N}
    a::A
end

to_complex_amplitude(V::TabulatedFieldPowerSpectrum, ϕ) =
    TabulatedFieldAmplitude(.√(V.a) .* exp.(im*ϕ))

@doc raw"""
    TabulatedVectorPotentialAmplitude(a)

Indicates that the values contained in the array `a` are of the kind
``A(q)``.
"""
struct TabulatedVectorPotentialAmplitude{T,N,A<:AbstractArray{T,N}} <: TabulatedFieldQuantity{T,N}
    a::A
end

function to_complex_amplitude(V::TabulatedVectorPotentialAmplitude, ϕ)
    all(iszero, ϕ) ||
    @warn "Non-zero phase for tabulated amplitudes meaningless, ignoring."
    V
end

@doc raw"""
    TabulatedVectorPotentialMagnitude(a)

Indicates that the values contained in the array `a` are of the kind
``|A(q)|``.
"""
struct TabulatedVectorPotentialMagnitude{T,N,A<:AbstractArray{T,N}} <: TabulatedFieldQuantity{T,N}
    a::A
end

to_complex_amplitude(V::TabulatedVectorPotentialMagnitude, ϕ) =
    TabulatedVectorPotentialAmplitude(V.a .* exp.(im*ϕ))

val_or_unity(x) = iszero(x) ? one(x) : x

function regrid(V, q; verbosity=0,
                qstart=minimum(q),
                spline_order=3, num_knots=min(100, length(q)), ξ=1.2, KnotSetType=LinearKnotSet,
                tol=√(eps(val_or_unity(maximum(abs, q)))), max_iter=30, kwargs...)
    Δqmin,Δqmax = extrema(abs, diff(q))
    rtol = tol*norm(V)

    verbosity > 0 &&
    @info "Regridding" V q Δqmin Δqmax spline_order num_knots ξ KnotSetType tol rtol max_iter

    qmin,qmax = extrema(q)
    qr = range(qstart, stop=qmax, step=Δqmin)
    @assert qmin < qmax # We cannot deal with infinitesimal intervals
    Bn(n) = BSpline(KnotSetType(spline_order, qmin, qmax, n))

    # p = sortperm(q)
    # qp = q[p]
    # Vp = selectdim(V, 1, p)

    Vre,Vim = real(V),imag(V)

    B = nothing
    C = nothing

    residuals = fill(Inf, max_iter)
    i = 0
    while i < max_iter
        i += 1
        B = Bn(num_knots)
        C = interpolate(B, q, Vre) + im*interpolate(B, q, Vim)
        residual = norm(V - (B[q,:]*C))
        verbosity > 2 && @info "Regridding, iteration $(i)/$(max_iter)" B residual
        residuals[i] = residual
        if i > 1 && residual > residuals[i-1] && verbosity > 1
            @warn "Residual increasing" residuals[1:i]
        end
        residual < rtol && break
        num_knots = ceil(Int, ξ*num_knots)
    end
    if residuals[i] > rtol
        @warn "Regridding did not succeed in $(max_iter) iterations"
    end

    B[qr,:]*C, qr
end

@doc raw"""
    time_resolved_vector_potential(F::TabulatedFieldAmplitude, t::AbstractVector{<:Unitful.Time})

Compute ``A\{t\}`` from ``F(t)`` (after potentially regridding to
uniform time grid) via `fft`, time-integrate in the Fourier domain
(multiply by ``-\frac{1}{\im\omega}``), and transform back using
`ifft`.
"""
function time_resolved_vector_potential(F::TabulatedFieldAmplitude, t::AbstractVector{<:Unitful.Time}; kwargs...)
    Fr,tr = regrid(F, t; kwargs...)
    time_resolved_vector_potential(TabulatedFieldAmplitude(Fr), tr; kwargs...)
end

function time_resolved_vector_potential(F::TabulatedFieldAmplitude, t::AbstractRange{<:Unitful.Time}; kwargs...)
    # Fourier transform
    @error "Not yet implemented"
    time_resolved_vector_potential(TabulatedFieldAmplitude(Fω), ωau; kwargs...)
end

# Perform time integral in Fourier domain to get vector potential from
# field amplitude.
function time_resolved_vector_potential(V::TabulatedFieldAmplitude, ωau::AbstractVector{<:AbstractFloat}; kwargs...)
    f = -inv.(im*ωau)
    A = f .* V.a
    time_resolved_vector_potential(TabulatedVectorPotentialAmplitude(A), ωau; kwargs...)
end

# Regrid to uniform ω grid
function time_resolved_vector_potential(A::TabulatedVectorPotentialAmplitude, ωau::AbstractVector; kwargs...)
    Ar,ωr = regrid(A, ωau; qstart=0, kwargs...)
    time_resolved_vector_potential(TabulatedVectorPotentialAmplitude(Ar), ωr; kwargs...)
end

to_angular_frequency_au_proportional(q, factor) =
    ustrip(q) * austrip(factor*unit(first(q)))

to_angular_frequency_au(ħω::AbstractVector{<:Unitful.Energy}) =
    to_angular_frequency_au_proportional(ħω, 1)

to_angular_frequency_au(f::AbstractVector{<:Unitful.Frequency}) =
    to_angular_frequency_au_proportional(f, 2π)

to_angular_frequency_au(f::AbstractVector{<:Unitful.Length}) =
    austrip.(2π*u"c"./f)

# Variable transform from q -> ω
time_resolved_vector_potential(V::Union{TabulatedFieldAmplitude,TabulatedVectorPotentialAmplitude},
                               q::AbstractVector{<:Quantity}; kwargs...) =
    time_resolved_vector_potential(V, to_angular_frequency_au(q); kwargs...)

# Finally, we can perform the inverse Fourier transform to get A{t}.
function time_resolved_vector_potential(A::TabulatedVectorPotentialAmplitude, ωau::AbstractRange; kwargs...)
    @assert all(≥(0), ωau)
    n = length(ωau)
    At = irfft(parent(A), 2n-1, 1)
    t = rfftfreq(length(ωau), 2π/(step(ωau)))
    TabulatedVectorPotentialAmplitude(At), t
end

time_resolved_vector_potential(A::TabulatedVectorPotentialAmplitude, t::AbstractVector{<:Unitful.Time}) = (A,t)

"""
    tabulated_field(q, V, ϕ=zeros(...); I₀=1)

Generate an electric field specified by the quantities tabulated: `q`
must be a `Unitful` vector (wavelengths, frequencies, energies, time,
&c.) with respect to which the [`TabulatedFieldQuantity`](@ref) `V` is
tabulated. The phase array `ϕ`, which must have the same dimensions as
`V`, default to a flat phase; it is assumed to be given in radians.

The resultant field will be normalized such that the integrated
intensity is `I₀` (cf. [Parseval's
theorem](https://en.wikipedia.org/wiki/Parseval's_theorem)).
"""
function tabulated_field(q::AbstractVector{<:Quantity},
                         V::TabulatedFieldQuantity{<:Any,N},
                         ϕ::AbstractArray{<:Any,N}=similar_fill(V, 0);
                         I₀=1, num_knots=100, spline_order=3,
                         kwargs...) where N
    Vz = to_complex_amplitude(V, ϕ)
    @info "Tabulated field" q V ϕ Vz

    # Vz is now either a TabulatedFieldAmplitude or a
    # TabulatedVectorPotentialAmplitude. In the former case, we have
    # to time-integrate (which we do in the Fourier domain), to get a
    # vector potential.

    A,t = time_resolved_vector_potential(Vz, q; kwargs...)

    # # Construct B-spline from A{t}
    # B = BSpline(LinearKnotSet(spline_order, austrip(first(t)), austrip(last(t)), num_knots))
    # verbosity > 1 && @info "Generated B-spline" num_knots B

    # BSplineField(B, t, A)
end

export TabulatedFieldAmplitude, TabulatedFieldMagnitude, TabulatedFieldPowerSpectrum,
    TabulatedVectorPotentialAmplitude, TabulatedVectorPotentialMagnitude,
    tabulated_field
