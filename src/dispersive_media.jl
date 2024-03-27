abstract type Medium <: DispersiveElement end

"""
    IsotropicMedium(material, d[, ∂k∂ω₀ = 0])

Describes the dispersion through an isotropic medium of thickness `d`
(such as an atomic gas), where the refractive index, given by
`material`, is the same in all directions.

The optional `∂k∂ω₀` can be used to subtract the linear slope from the
dispersion relation, i.e. remove the time shift, such that a pulse
with central angular frequency `ω₀` stays centred in the frame of
reference.
"""
struct IsotropicMedium{Material,D<:Length,U} <: Medium
    material::Material
    d::D
    ∂k∂ω₀::U
end
IsotropicMedium(material, d; ω₀=nothing) =
    IsotropicMedium(material, d, isnothing(ω₀) ? 0 : dispersion_slope(material, ω₀))

tropy(::IsotropicMedium) = Isotropic()

function Base.show(io::IO, m::IsotropicMedium)
    dau = austrip(m.d)
    du = u"μm"
    printfmt(io, "IsotropicMedium({1:0.2g} Bohr = {2:0.2g} {3:s} of ",
             dau, ustrip(auconvert(du, dau)), du)
    show(io, m.material)
    printfmt(io, ", ∂k∂ω₀ = {1:0.4g})", m.∂k∂ω₀)
end

function frequency_response(m::IsotropicMedium, ω::Number)
    f = ω/2π
    n = real(m.material(auconvert(u"Hz", f)))

    # Wavevector in vacuum
    k₀ = ω/austrip(1u"c")
    # Wavevector in the medium, subtracting the linear component to
    # centre the pulse in the frame of reference.
    k = n*k₀ - ω*m.∂k∂ω₀

    d = austrip(m.d)

    exp(-im*k*d)
end

"""
    Crystal(material, d, R[, ∂k∂ω₀ = 0])

Describes the dispersion through a general crystal of thickness `d`
and orientation given by the rotation `R`, with different refractive
indices, given by the components of `material`, along the different
crystal axes.

In the special case where `material` has two components, we have a
_uniaxial_ crystal, which is birefringent, and `material[1]` is
referred to as _ordinary_ the refractive index, and `material[2]` as
the _extraordinary_ refractive index.

The optional `∂k∂ω₀` can be used to subtract the linear slope from the
dispersion relation, i.e. remove the time shift, such that a pulse
with central angular frequency `ω₀` stays centred in the frame of
reference.
"""
struct Crystal{Material,D<:Length,Rotation,U} <: Medium
    material::Material
    d::D
    R::Rotation
    ∂k∂ω₀::U
end
Crystal(material, d, R=I; ω₀=nothing) =
    Crystal(material, d, R,
            isnothing(ω₀) ? 0 : dispersion_slope(material, ω₀))

function Base.show(io::IO, m::Crystal)
    dau = austrip(m.d)
    du = u"μm"
    printfmt(io, "Crystal({1:0.2g} Bohr = {2:0.2g} {3:s} of ",
             dau, ustrip(auconvert(du, dau)), du)
    show(io, m.material)
    printfmt(io, ", R = {1:s}, ∂k∂ω₀ = {2:0.4g})", m.R, m.∂k∂ω₀)
end

function frequency_response(m::Crystal, ω::Number)
    SVector(1,1,1)
end

export IsotropicMedium, Crystal
