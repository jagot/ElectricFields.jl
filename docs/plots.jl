using ElectricFields
using Unitful

using PyPlot
using Jagot.plotting
plot_style("ggplot")

function savedocfig(name,dir="figures")
    fig = gcf()
    filename = joinpath(@__DIR__, "src", dir, "$(name).svg")
    savefig(filename,
            transparent=false,
            facecolor=fig.get_facecolor())
    close(fig)
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

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
mkpath("docs/src/figures")
@echo index_example()
@echo index_polarized_example()
