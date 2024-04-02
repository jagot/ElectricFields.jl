import ElectricFields: CascadedDispersiveElement, frequency_response

@testset "Dispersed fields" begin
    @testset "Chirped fields" begin
        function chirped_gaussian(t::Number, ω₀, ϕ₀, τ, η; verbose = false)
            γ = τ^2/(8log(2))
            g = inv(γ^2 + η^2)

            ϕη = atan(-η,γ)/2
            A = √(γ*√g)*exp(im*ϕη)

            aη = γ*g/4
            bη = η*g/2

            verbose && @info "Pulse params" ω₀ γ η ϕ₀ ϕη A aη bη

            ca = exp(im*(ω₀*t + ϕ₀ + bη/2*t^2))
            env = exp(-aη*t^2)

            ca′ = im*(ω₀ + bη*t)*exp(im*(ω₀*t + ϕ₀ + bη/2*t^2))
            env′ = (-2aη*t)*exp(-aη*t^2)
            1/ω₀*imag(A*ca*env), -1/ω₀*imag(A*(ca′*env + ca*env′))
        end

        function chirped_gaussian(t::AbstractVector, args...; kwargs...)
            T = eltype(t)
            if !isempty(t)
                A,F = chirped_gaussian(first(t), args...; kwargs...)
                T = promote_type(T, typeof(A), typeof(F))
            end
            V = zeros(T, length(t), 2)
            for (i,t) in enumerate(t)
                V[i,:] .= chirped_gaussian(t, args...; kwargs...)
            end
            V[:,1],V[:,2]
        end

        @testset "$(name) pulse" for (name,τ,ηs,rtol) in (
            ("Short", 3u"fs", 5.0u"fs^2"*[1,-1], 4e-2),
            ("Long", 30u"fs", 15.0u"fs^2"*[1,-1], 5e-3))
            @testset "η = $(η)" for η in ηs
                τ = austrip(τ)
                η = austrip(η)

                @field(F) do
                    λ = 800u"nm"
                    I₀ = 1.0
                    τ = τ
                    σoff = 4.0
                    σmax = 6.0
                    env = :trunc_gauss
                end

                ω₀ = photon_energy(F)
                Fc = chirp(F, η, ω₀)

                torig = timeaxis(F)
                t = timeaxis(Fc)

                @test t[1] ≲ torig[1] rtol=5e-3
                @test t[end] ≳ torig[end] rtol=5e-3

                Fv = field_amplitude(Fc, t)
                Av = vector_potential(Fc, t)

                F₀ = amplitude(F)
                Avref,Fvref = F₀ .* chirped_gaussian(t, ω₀, 0.0, τ, η)

                test_approx_eq(Fv, Fvref, rtol=rtol)
                test_approx_eq(Av, Avref, rtol=rtol)
            end
        end
    end

    @testset "Cascaded dispersive elements" begin
        a = PhaseShift(π)
        b = Chirp(austrip(5u"fs^2"), 1.0)
        c = CascadedDispersiveElement(())

        d = a*b
        e = d*b
        f = a*d
        g = d*d

        @test frequency_response(c, 0.4) == 1

        @test d == CascadedDispersiveElement((a,b))
        @test frequency_response(d, 0.4) ≈ frequency_response(a, 0.4)*frequency_response(b, 0.4)

        @test e == CascadedDispersiveElement((a,b,b))
        @test frequency_response(e, 0.4) ≈ frequency_response(a, 0.4)*frequency_response(b, 0.4)^2

        @test f == CascadedDispersiveElement((a,a,b))
        @test frequency_response(f, 0.4) ≈ frequency_response(a, 0.4)^2*frequency_response(b, 0.4)

        @test g == CascadedDispersiveElement((a,b,a,b))
        @test frequency_response(g, 0.4) ≈ frequency_response(a, 0.4)^2*frequency_response(b, 0.4)^2
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
        Fc2 = phase_shift(Fc, π)

        torig = timeaxis(F)
        t = timeaxis(Fc)

        @test t[1] ≲ torig[1] rtol=5e-3
        @test t[end] ≳ torig[end] rtol=5e-3

        @test field_amplitude(Fc2, t) ≈ -field_amplitude(Fc, t)

        @test string(Fc) == "DispersedField($(F), $(de))"
        @test string(Fc2) == "DispersedField($(F), $(Fc2.de))"

        @test pretty_print_object(Fc) == "DispersedField:\n$(pretty_print_object(F))\n  – dispersed through $(pretty_print_object(Fc.de))"
        @test pretty_print_object(Fc2) == "DispersedField:\n$(pretty_print_object(F))\n  – dispersed through $(pretty_print_object(Fc2.de))"
    end

    @testset "Crystal($(material))" for material in (BK7,Quartz,KTP)
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
        de = Crystal(material, 12u"μm", ω₀=ω₀)
        Fc = DispersedField(F, de, verbosity=4)
        Fc2 = phase_shift(Fc, π)

        torig = timeaxis(F)
        t = timeaxis(Fc)

        @test t[1] ≲ torig[1] rtol=5e-3
        @test t[end] ≳ torig[end] rtol=5e-3

        @test field_amplitude(Fc2, t) ≈ -field_amplitude(Fc, t)

        @test string(Fc) == "DispersedField($(F), $(de))"
        @test string(Fc2) == "DispersedField($(F), $(Fc2.de))"

        @test pretty_print_object(Fc) == "DispersedField:\n$(pretty_print_object(F))\n  – dispersed through $(pretty_print_object(Fc.de))"
        @test pretty_print_object(Fc2) == "DispersedField:\n$(pretty_print_object(F))\n  – dispersed through $(pretty_print_object(Fc2.de))"
    end
end
