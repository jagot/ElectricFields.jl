@testset "Field spectra" begin
    @testset "Gaussian, $(kind)" for kind = [:linear, :transverse, :elliptical]
        if kind == :elliptical
            @field(F) do
                I₀ = 1.0
                T = 1.0
                τ = 2.0
                σmax = 6.0
                ξ = 0.67
            end
        else
            @field(F) do
                I₀ = 1.0
                T = 1.0
                τ = 2.0
                σmax = 6.0
                kind = kind
            end
        end

        t = timeaxis(F)
        ω = fftshift(fftω(t))

        # We need to undo the phase, since the FFT does not care that
        # pulse is centred around zero.
        F̂v = exp.(im*ω*t[1]) .* fftshift(nfft(F, t), 1)
        Âv = exp.(im*ω*t[1]) .* fftshift(nfft_vector_potential(F, t), 1)

        sel = kind == :linear ? 1 : (kind == :transverse ? 3 : Colon())
        @test norm(F̂v[:,sel]) ≈ 1.4251872372067347
        @test norm(Âv[:,sel]) ≈ 0.22581870287199948

        F̂v_exact = field_amplitude_spectrum(F, ω)
        Âv_exact = vector_potential_spectrum(F, ω)

        test_approx_eq(F̂v, F̂v_exact, rtol=0.2)
        test_approx_eq(abs2.(F̂v), abs2.(F̂v_exact), rtol=1e-2)

        test_approx_eq(Âv, Âv_exact, rtol=0.2)
        test_approx_eq(abs2.(Âv), abs2.(Âv_exact), rtol=1e-2)

        test_approx_eq(field_amplitude_spectrum(delay(F, 6.0), ω), exp.(-im*6ω) .* F̂v_exact)
    end

    @testset "Constant field" begin
        tmax = 5.0

        @field(F) do
            E₀ = 1.0
            tmax = tmax
            kind = :constant
        end

        t = range(0, stop=20tmax, length=2001)
        ω = fftshift(fftω(t))

        # We need to undo the phase, since the FFT does not care that
        # pulse is centred around zero.
        F̂v = exp.(im*ω*t[1]) .* fftshift(nfft(F, t))
        F̂v_exact = field_amplitude_spectrum(F, ω)

        # The vector potential seems broken
        Âv = exp.(im*ω*t[1]) .* fftshift(nfft_vector_potential(F, t))
        Âv_exact = vector_potential_spectrum(F, ω)

        # We do not yet have phase agreement
        test_approx_eq(F̂v, F̂v_exact, rtol=0.2, isbroken=true)
        test_approx_eq(abs2.(F̂v), abs2.(F̂v_exact), rtol=1e-2)

        test_approx_eq(Âv, Âv_exact, rtol=0.2, isbroken=true)
        test_approx_eq(abs2.(Âv), abs2.(Âv_exact), rtol=1e-2, isbroken=true)
    end
end
