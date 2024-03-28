@testset "Dispersed fields" begin
    @testset "Chirped fields" begin
        # This reference function produces chirped _electric fields_,
        # whereas ElectricFields.jl compute the electric field from
        # the _vector potential_. This will impact the comparison for
        # pulses of short durations.
        function chirped_gaussian(t, ω₀, ϕ₀, τ, η; verbose = false)
            γ = τ^2/(8log(2))
            g = inv(γ^2 + η^2)

            ϕη = atan(-η,γ)/2
            A = √(γ*√g)*exp(im*ϕη)

            aη = γ*g/4
            bη = η*g/2

            verbose && @info "Pulse params" ω₀ γ η ϕ₀ ϕη A aη bη

            ca = exp(im*(ω₀*t + ϕ₀ + bη/2*t^2))
            env = exp(-aη*t^2)
            real(A*ca*env)
        end

        @testset "Short pulse, η = $(η)" for η = 5.0u"fs^2"*[1,-1]
            τ = austrip(3u"fs")
            η = austrip(η)

            @field(F) do
                λ = 800u"nm"
                I₀ = 1.0
                τ = τ
                σoff = 4.0
                σmax = 6.0
                env = :trunc_gauss
                ϕ = π
            end

            ω₀ = photon_energy(F)
            Fc = chirp(F, η, ω₀)

            torig = timeaxis(F)
            t = timeaxis(Fc)

            @test t[1] < torig[1]
            @test t[end] > torig[end]

            Fv = field_amplitude(Fc, t)
            Fvref = chirped_gaussian.(t, ω₀, 0.0, τ, η)

            test_approx_eq(Fv, Fvref, rtol=3e-1)
        end

        @testset "Long pulse, η = $(η)" for η = 15.0u"fs^2"*[1,-1]
            τ = austrip(30u"fs")
            η = austrip(η)

            @field(F) do
                λ = 800u"nm"
                I₀ = 1.0
                τ = τ
                σoff = 4.0
                σmax = 6.0
                env = :trunc_gauss
                ϕ = π
            end

            ω₀ = photon_energy(F)
            Fc = chirp(F, η, ω₀)

            torig = timeaxis(F)
            t = timeaxis(Fc)

            @test t[1] ≲ torig[1] rtol=5e-3
            @test t[end] ≳ torig[end] rtol=5e-3

            Fv = field_amplitude(Fc, t)
            Fvref = chirped_gaussian.(t, ω₀, 0.0, τ, η)

            test_approx_eq(Fv, Fvref, rtol=3e-2)
        end
    end

    @testset "Isotropic material" begin
        τ = austrip(3u"fs")

        @field(F) do
            λ = 800u"nm"
            I₀ = 1.0
            τ = τ
            σoff = 4.0
            σmax = 12.0
            env = :trunc_gauss
            ϕ = π
        end

        ω₀ = photon_energy(F)
        de = IsotropicMedium(BK7, 220u"μm", ω₀=ω₀)
        Fc = DispersedField(F, de, verbosity=4)

        torig = timeaxis(F)
        t = timeaxis(Fc)

        @test t[1] ≲ torig[1] rtol=5e-3
        @test t[end] ≳ torig[end] rtol=5e-3
    end
end
