"""
    span(f)

Return the end-points for the [`envelope`](@ref) of the field `f`.
"""
span(f::AbstractField) = span(envelope(f))

const DEFAULT_SAMPLING_FACTOR = 100
default_sampling_frequency(f::AbstractField) = DEFAULT_SAMPLING_FACTOR*austrip(max_frequency(f))

@doc raw"""
    steps(f[, fs])

Return number of time steps for the field `f`, optionally supplying a
sampling frequency.

We construct the time axis such that the highest frequency component
of the field is resolved. By the Nyquist sampling theorem, we need
minimum ``f_s>2f_{\textrm{max}}``, but to be on the safe side, we use,
as default, ``f_s=100f_{\textrm{max}}``. This also makes plots nicer.

If `fs<:Integer`, the sampling frequency is computed as `fs/T`,
i.e. `fs` steps per cycle.
"""
steps(f::AbstractField, fs::Real=default_sampling_frequency(f); s=span(f)) =
    ceil(Int, fs*abs(-(endpoints(s)...)))
steps(f::AbstractField, ndt::Int; kwargs...) = steps(f, ndt/austrip(period(f)); kwargs...)

"""
    timeaxis(f[, fs])

Construct a time axis for the field `f`, covering the [`span`](@ref)
of the envelope in the number of [`steps`](@ref) given by the sample
frequency `fs`.
"""
function timeaxis(f::AbstractField, fs::Number=default_sampling_frequency(f); s=span(f))
    num_steps = steps(f, fs; s=s)
    if num_steps > 1 || s.left == s.right
        range(s, length=num_steps)
    else
        # This hack is mainly useful for calculations where exactly
        # one time step should be taken, (b-a) becomes the step
        # length.
        a,b = endpoints(s)
        a:(b-a):a
    end
end


"""
    timeaxis(f, (N,dt))

Alternate interface that constructs the time axis with `N` steps
starting from `first(span(f))`, with intervals of `dt`.
"""
function timeaxis(f::AbstractField, (N,dt)::Tuple{<:Integer,<:Real})
    a,b = endpoints(span(f))
    t = range(a, length=N, step=dt)
    t[end] < b && @warn "$(N) steps of $(dt) jiffies does not cover the span of the field, $(b) ∉ $(t)"
    t
end

export span, steps, timeaxis
