"""
    BSplineField(B, C)

Represents an electric field using an expansion over a B-spline basis
`B`, with `C` being the array (vector for 1d fields, 3-column matrix
for 3d fields) expansion coefficients for the vector potential.
"""
struct BSplineField{Bt,Ct<:AbstractArray} <: AbstractField
    B::Bt
    C::Ct
end

function BSplineField(B::BSplineOrView, t::AbstractVector, A::AbstractVecOrMat)
    BSplineField(B, interpolate(B, t, A))
end

BSplineField(B::BSplineOrView, t::AbstractVector{<:Unitful.Time}, A::AbstractVecOrMat) =
    BSplineField(B, austrip.(t), A)

Base.show(io::IO, f::BSplineField) =
    write(io, "B-spline field")

function Base.show(io::IO, ::MIME"text/plain", f::BSplineField)
    show(io, f)
    println(io, " expanded over")
    write(io, "  ")
    show(io, f.B)
end

dimensions(::BSplineField{<:Any,<:AbstractVector}) = 1
dimensions(::BSplineField{<:Any,<:AbstractMatrix}) = 3

polarization(::BSplineField{<:Any,<:AbstractVector}) = LinearPolarization()
polarization(::BSplineField{<:Any,<:AbstractMatrix}) = ArbitraryPolarization()

linear_combination(V::AbstractVector, C::AbstractVector) =
    dot(V, C)

linear_combination(V::AbstractVector, C::AbstractMatrix) =
    SVector(dot(V, view(C, :, 1)), dot(V, view(C, :, 2)), dot(V, view(C, :, 3)))

vector_potential(f::BSplineField, t::Number) =
    linear_combination(f.B[t,:], f.C)

vector_potential(f::BSplineField, t::AbstractVector) =
    f.B[t,:]*f.C

field_amplitude(f::BSplineField, t::Number) =
    -linear_combination(derivative(f.B, t, :, 1), f.C)

field_amplitude(f::BSplineField, t::AbstractVector) =
    -derivative(f.B, t, :, 1)*f.C

export BSplineField
