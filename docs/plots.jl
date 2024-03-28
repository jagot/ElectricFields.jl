using ElectricFields
using Unitful
using UnitfulAtomic
using FFTW
using Statistics

using PythonPlot
using Jagot
using Jagot.plotting
plot_style("ggplot")

function wind_transf(t, x, w)
    f = rfftfreq(length(t),1.0/(t[2]-t[1]))
    ecat(a,b) = cat(a,b,dims=ndims(x)+1)
    Cₓ = reduce(ecat, x.*circshift(w,i) for i = 1:length(t))
    Wₓ = rfft(Cₓ, 1)
    f,Wₓ
end
gabor(t, x, σ) = wind_transf(t, x, fftshift(exp.(-(t .- mean(t)).^2/2σ^2)))

function savedocfig(name,dir="figures")
    fig = gcf()
    filename = joinpath(@__DIR__, "src", dir, "$(name).svg")
    savefig(filename,
            transparent=false,
            facecolor=fig.get_facecolor())
    PythonPlot.close("all")
    if isfile(filename)
        println("Saved $(name) to $(filename)")
    else
        @warn "Saving $(name) to $(filename) failed"
    end
end

function savedfigure(fun::Function, name, args...; kwargs...)
    cfigure(name, args...; kwargs...) do
        fun()
    end
    savedocfig(name)
end

function index_example()
    @field(IR) do
        I₀ = 1e14u"W/cm^2"
        λ = 800.0u"nm"
        τ = 6.2u"fs"
        σmax = 6.0
    end
    t = timeaxis(IR)
    tplot = 24.3e-3t
    Fv = field_amplitude(IR, t)
    Fenv = field_envelope.(IR, t)
    Av = vector_potential(IR, t)
    Iv = instantaneous_intensity.(IR, t)
    Ienv = intensity.(IR, t)

    savedfigure("index_example") do
        csubplot(311, nox=true) do
            plot(tplot, Av)
            ylabel(L"$A(t)$ [au]")
        end
        csubplot(312, nox=true) do
            plot(tplot, Fv)
            plot(tplot, Fenv)
            ylabel(L"$F(t)$ [au]")
        end
        csubplot(313) do
            plot(tplot, Iv)
            plot(tplot, Ienv)
            ylabel(L"$I(t)$ [au]")
        end
        xlabel(L"$t$ [fs]")
    end
end

function index_polarized_example()
    @field(A) do
        I₀ = 1.0
        ω = 1.0
        cycles = 6.0
        env = :cos²
        ξ = 1.0
    end

    @field(B) do
        I₀ = 1.0
        ω = 2.0
        cycles = 6.0
        env = :cos²
        ξ = -1.0
    end

    F = A + delay(B, 3/2π)

    t = timeaxis(F)
    Fv = field_amplitude(F, t)

    savedfigure("index_polarized_example", figsize=(8,10)) do
        csubplot(211) do
            plot(t, Fv[:,1], label=L"F_x(t)")
            plot(t, Fv[:,3], label=L"F_z(t)")
            legend()
            xlabel(L"$t$ [jiffies]")
        end
        csubplot(212, projection="3d") do
            plot3D(t, Fv[:,1], Fv[:,3])
            xlabel(L"$t$ [jiffies]")
            ylabel(L"F_x(t)")
            gca().set_zlabel(L"F_z(t)")
        end
    end
end

function index_spectrum_example()
    @field(F) do
        I₀ = 1.0
        T = 1.0
        τ = 2.0
        σmax = 6.0
    end

    t = timeaxis(F)
    Fv = field_amplitude(F, t);
    Av = vector_potential(F, t);

    ω = fftshift(fftω(t));
    # We need to undo the phase, since the FFT does not care that
    # pulse is centred around zero.
    F̂v = exp.(im*ω*t[1]) .* fftshift(nfft(F, t), 1);
    Âv = exp.(im*ω*t[1]) .* fftshift(nfft_vector_potential(F, t), 1);

    F̂v_exact = field_amplitude_spectrum(F, ω);
    Âv_exact = vector_potential_spectrum(F, ω);

    sel = ind(ω, -20):ind(ω, 20)

    savedfigure("index_spectrum_example", figsize=(8,10)) do
        csubplot(311) do
            plot(t, Fv, label=L"F(t)")
            plot(t, Av, label=L"A(t)")
            axes_labels_opposite(:x)
            xlabel(L"t/T")
            legend(loc=1)
        end
        csubplot(323, nox=true) do
            semilogy(ω[sel], abs2.(F̂v[sel,:]), label=L"$|F(\omega)|^2$, FFT")
            semilogy(ω[sel], abs2.(Âv[sel,:]), label=L"$|A(\omega)|^2$, FFT")
            ax = axis()
            semilogy(ω[sel], abs2.(F̂v_exact[sel,:]), "--", label=L"$|F(\omega)|^2$, exact")
            semilogy(ω[sel], abs2.(Âv_exact[sel,:]), "--", label=L"$|A(\omega)|^2$, exact")
            axis(ax)
            legend(loc=3)
        end
        csubplot(325) do
            semilogy(ω[sel], abs.(abs2.(F̂v[sel,:])-abs2.(F̂v_exact[sel,:])), label=L"||\hat{F}_{\mathrm{FFT}}|^2-|\hat{F}_{\mathrm{exact}}|^2|")
            semilogy(ω[sel], abs.(abs2.(Âv[sel,:])-abs2.(Âv_exact[sel,:])), label=L"||\hat{A}_{\mathrm{FFT}}|^2-|\hat{A}_{\mathrm{exact}}|^2|")
            xlabel(L"$\omega$ [rad/jiffies]")
            legend(loc=3)
        end
        csubplot(324, nox=true) do
            plot(ω[sel], (angle.(F̂v[sel,:])), label=L"$\arg\{F(\omega)\}$, FFT")
            plot(ω[sel], (angle.(Âv[sel,:])), label=L"$\arg\{A(\omega)\}$, FFT")
            plot(ω[sel], (angle.(F̂v_exact[sel,:])), "--", label=L"$\arg\{F(\omega)\}$, exact")
            plot(ω[sel], (angle.(Âv_exact[sel,:])), "--", label=L"$\arg\{A(\omega)\}$, exact")
            axes_labels_opposite(:y)
            legend(loc=3)
            π_labels(:y)
        end
        csubplot(326) do
            semilogy(ω[sel], abs2.(F̂v[sel,:] - F̂v_exact[sel,:]), label=L"|\hat{F}_{\mathrm{FFT}}-\hat{F}_{\mathrm{exact}}|")
            semilogy(ω[sel], abs2.(Âv[sel,:] - Âv_exact[sel,:]), label=L"|\hat{A}_{\mathrm{FFT}}-\hat{A}_{\mathrm{exact}}|")
            axes_labels_opposite(:y)
            xlabel(L"$\omega$ [rad/jiffies]")
            legend(loc=3)
        end
    end
end

function apodized_field()
    @field(F) do
        ω = 1.0
        I₀ = 1.0
        ramp = 0.0
        flat = 3.0
        env = :tophat
    end

    t = timeaxis(F)
    tplot = 24.2e-3t

    tmin = 1.0
    tmax = 14.0

    Fw = ApodizedField(F, tmin, tmax, ElectricFields.Kaiser(3))

    Fv = field_amplitude(F, t)
    Av = vector_potential(F, t)

    Fwv = field_amplitude(Fw, t)
    Awv = vector_potential(Fw, t)

    w = ElectricFields.window_value.(Fw.window, 1.0, 14.0, t)

    savedfigure("apodized_field", figsize=(7,8)) do
        csubplot(211, nox=true) do
            plot(tplot, Fv)
            plot(tplot, Fwv)
            ylabel(L"F(t)")
        end
        csubplot(212) do
            plot(tplot, Av, label="Original field")
            plot(tplot, Awv, label="Windowed field")
            plot(tplot, w, "--", label="Window")
            legend(loc=1)
            xlabel(L"$t$ [fs]")
            ylabel(L"A(t)")
        end
    end
end

function apodizing_windows()
    x = range(-0.55, stop=0.55, length=1001)

    savedfigure("windows", figsize=(7,8)) do
        ws = (ElectricFields.Hann(), ElectricFields.Hamming(),
              ElectricFields.Blackman(),
              ElectricFields.BlackmanExact(),
              ElectricFields.Nuttall(),
              ElectricFields.BlackmanNuttall(),
              ElectricFields.BlackmanHarris(),
              ElectricFields.Kaiser(3), ElectricFields.Kaiser(2),
              ElectricFields.Rect())

        csubplot(211, nox=true) do
            for w in ws
                plot(x, ElectricFields.window_value.(w, x), label=string(w))
            end
            ylabel(L"w(x)")
        end
        csubplot(212) do
            for w in ws
                plot(x, ElectricFields.window_derivative.(w, x), label=string(w))
            end
            legend(loc=1, ncol=2)
            xlabel(L"$t$ [fs]")
            ylabel(L"w'(x)")
        end
    end
end

function chirped_field()
    τ = austrip(3u"fs")
    η = austrip(5.0u"fs^2")

    @field(F) do
        λ = 800u"nm"
        I₀ = 1.0
        τ = τ
        σoff = 4.0
        σmax = 6.0
        env = :trunc_gauss
        ϕ = π
    end
    Fc = chirp(F, η)

    γ = τ^2/(8log(2))
    ω₀ = photon_energy(F)

    t = range(span(Fc), length=1500)
    tplot = ustrip(auconvert.(u"fs", 1))*t

    Av = vector_potential(F, t)
    Fv = field_amplitude(F, t)

    Frec = @time field_amplitude(Fc, t)
    Arec = @time vector_potential(Fc, t)

    freq,G = gabor(t, Frec, 1austrip(period(F)))
    fsel = ind(freq, 0):ind(freq, 2ω₀/2π)

    savedfigure("chirped_field", figsize=(7,8)) do
        csubplot(311) do
            plot(tplot, Av, "k")
            plot(tplot, Arec)
            xlabel(L"$t$ [fs]")
            axes_labels_opposite(:x)
            ylabel(L"A(t)")
        end
        csubplot(312,nox=true) do
            plot(tplot, Fv, "k")
            plot(tplot, Frec)
            ylabel(L"F(t)")
        end
        csubplot(313) do
            plot_map(tplot, freq[fsel], abs2.(G[fsel,:]))
            hline(ω₀/2π, color="white")
            yl = ylim()
            plot(tplot, (ω₀ .+ 1/2*η/(γ^2 + η^2)*t)/2π, color="white")
            ylim(yl)
            xlabel(L"$t$ [fs]")
            ylabel(L"\omega/2\pi")
        end
    end
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
fig_dir = joinpath(@__DIR__, "src", "figures")
mkpath(fig_dir)
@echo index_example()
@echo index_polarized_example()
@echo index_spectrum_example()
@echo apodized_field()
@echo apodizing_windows()
@echo chirped_field()
