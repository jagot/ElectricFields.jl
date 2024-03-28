"""
    DispersedField(f, de, B, s)

Represents the field `f` dispersed through the
[`DispersiveElement`](@ref) `de`. Since the dispersion is most
efficiently calculated in the frequency domain, this is done via
[`rfft`](@ref)/[`irfft`](@ref), and the dispersed field is
interpolated as the [`BSplineField`](@ref) `B`. The temporal [`span`](@ref)
of the dispersed field is given by `s`.
"""
struct DispersedField{Ft,Dt,Bt,S} <: AbstractField
    f::Ft
    de::Dt
    FB::Bt
    s::S
end

"""
    DispersedField(f, de; spline_order=7, Bfs=20/period(F))

Convenience constructor for [`DispersedField`](@ref) that
automatically finds the appropriate time span for `f` after it has
been dispersed through `de`, and then fits a [`BSpline`](@ref) to the
resultant vector potential, such that [`vector_potential`](@ref) and
[`field_amplitude`](@ref) can be evaluated at arbitrary times. The
number of knots can be influenced by the \"sampling frequency\" `Bfs`.
"""
function DispersedField(f::AbstractField, de;
                        spline_order=3, Bfs=20/period(parent(f)),
                        cutoff=1e5*√(eps()),
                        verbosity=0, kwargs...)
    s = find_time_span(f, de; verbosity=verbosity-2, kwargs...)
    t = timeaxis(f, s=s)
    s₀ = span(f)

    ω = rfftω(t)
    H = frequency_response(de, ω)

    As = rfft_vector_potential(f, t)
    Arec = irfft(H.*As, t)

    # Truncate time span such that |A| > cutoff*max|A|, but still
    # cover the whole time span of the undispersed field.
    if cutoff > 0
        Amag = vec(sum(abs, Arec, dims=2))
        Amax = maximum(Amag)
        abs_cutoff = cutoff*Amax

        sel = findall(∈(s₀), t)
        a = min(findfirst(>(abs_cutoff), Amag), first(sel))
        b = max(findlast(>(abs_cutoff), Amag), last(sel))

        s = t[a] .. t[b]
        t = t[a:b]
        Arec = selectdim(Arec, 1, a:b)

        verbosity > 1 && @info "Truncated to time interval $(s)" a b cutoff abs_cutoff
    end

    num_knots = ceil(Int, austrip(Bfs)*(s.right-s.left))
    B = BSpline(LinearKnotSet(spline_order, s.left, s.right, num_knots))
    verbosity > 1 && @info "Generated B-spline" num_knots B

    FB = BSplineField(B, t, Arec)

    DispersedField(f, de, FB, s)
end

@doc raw"""
    DispersedField(f::DispersedField, de)

Stacking another [`DispersiveElement`](@ref) `de` onto an already
dispersed field `f` is accomlished by multiplying the transfer
functions together, with the earliest applied as the right-most factor:
```math
H(\omega) = H_n(\omega) ... H_3(\omega)H_2(\omega)H_1(\omega).
```
"""
DispersedField(f::DispersedField, de; kwargs...) =
    DispersedField(parent(f), de*f.de; kwargs...)

function show(io::IO, f::DispersedField)
    write(io, "DispersedField(")
    show(io, f.f)
    write(io, ", ")
    show(io, f.de)
    write(io, ")")
end
function show(io::IO, mime::MIME"text/plain", f::DispersedField)
    println(io, "DispersedField:")
    show(io, mime, f.f)
    write(io, "\n  – dispersed through ")
    show(io, mime, f.de)
end

Base.parent(f::DispersedField) = f.f

for fun in [:vector_potential, :field_amplitude]
    @eval $fun(f::DispersedField, t::Number) =
        $fun(f.FB, t)
end

# Some of these forwards are iffy, i.e. after e.g. a chirp, the
# carrier is most certainly not a FixedCarrier anymore, and the
# duration is not the Fourier-limited one.
for fun in [:params, :carrier, :envelope, :polarization,
            :wavelength, :period, :frequency, :max_frequency,
            :wavenumber, :fundamental, :photon_energy,
            :intensity, :amplitude, :duration, :continuity,
            :dimensions, :rotation_matrix]
    @eval $fun(f::DispersedField) = $fun(parent(f))
end

span(f::DispersedField) = f.s

export DispersedField
