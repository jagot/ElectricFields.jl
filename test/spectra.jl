@testset "Field spectra" begin
    @testset "Gaussian" begin
        @field(F) do
            I₀ = 1.0
            T = 1.0
            τ = 2.0
            σmax = 6.0
        end

        t = timeaxis(F)
        Fv = field_amplitude(F, t)

        ω = fftshift(fftω(t))
        # We need to undo the phase, since the FFT does not care that
        # pulse is centred around zero.
        F̂v = -exp.(im*ω*t[1]) .* fftshift(nfft(F, t))
        F̂v_exact = field_amplitude_spectrum(F, ω)

        test_approx_eq(F̂v, F̂v_exact, rtol=0.2)
        test_approx_eq(abs2.(F̂v), abs2.(F̂v_exact), rtol=1e-2)
    end
end
