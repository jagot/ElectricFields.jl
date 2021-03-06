#   We construct the time axis such that the highest frequency component
#   of the field is resolved. By the Nyquist sampling theorem, we need
#   minimum \(f_s>2f_{\textrm{max}}\), but to be on the safe side, we
#   use, as default, \(f_s=100f_{\textrm{max}}\). This also makes plots
#   nicer.

span(f::AbstractField) = span(envelope(f))

const DEFAULT_SAMPLING_FACTOR = 100
default_sampling_frequency(f::AbstractField) = DEFAULT_SAMPLING_FACTOR*max_frequency(f)

steps(f::AbstractField, fs::Unitful.Frequency=default_sampling_frequency(f)) =
    ceil(Int, fs*abs(-(span(f)...)))
steps(f::AbstractField, ndt::Integer) = steps(f, ndt/period(f))

function timeaxis(f::AbstractField, fs::Unitful.Frequency=default_sampling_frequency(f))
    a,b = span(f)
    num_steps = steps(f, fs)
    if num_steps > 1
        range(a, stop=b, length=num_steps)
    else
        # This hack is mainly useful for calculations where exactly
        # one time step should be taken, (b-a) becomes the step
        # length.
        a:(b-a):a
    end
end
timeaxis(f::AbstractField, ndt::Integer) = timeaxis(f, ndt/period(f))

export span, steps, timeaxis
