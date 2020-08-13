@testset "Sellmeier" begin
    @testset "String representation" begin

        Bstr = string(BK7.B)
        Cstr = join(string.(BK7.C), ", ")
        @test string(BK7) == "Medium($(Bstr), [$(Cstr)])"# 
    end

    @test BK7(0.5876u"μm") ≈ 1.5168 atol=1e-3
    @test SiO2(0.5876u"μm") ≈ 1.4585 atol=1e-3# 

    @testset "Dispersion" begin
        # Dispersing a pulse through a positive amount of glass should
        # /delay/ the pulse, i.e. its maximum should arrive /later/,
        # and vice versa for a negative amount of glass (achievable
        # through precompensation, common in pulse compressors).

        λ = 500.0u"nm"
        f₀ = u"c"/λ |> u"THz"
        ω₀ = 2π*f₀
        τ = 6.2u"fs" # Pulse duration, intensity FWHM
        γ = τ^2/8log(2)

        f = range(0,stop=30,length=2000)*f₀
        ω = 2π*f

        Ê = exp.(-(ω .- ω₀).^2*γ)
        Ê′ = Ê.*dispersion(BK7, 6u"μm", f)
        Ê′′ = Ê.*dispersion(BK7, -6u"μm", f)
        Ê′′′ = Ê.*dispersion(BK7, -6u"μm", f, f₀)

        time_domain_envelope(spectrum) = abs.(fftshift(ifft(spectrum)*√(length(spectrum))))

        @test argmax(time_domain_envelope(Ê′)) > argmax(time_domain_envelope(Ê))
        @test argmax(time_domain_envelope(Ê′′)) < argmax(time_domain_envelope(Ê))
        @test argmax(time_domain_envelope(Ê′′′)) == argmax(time_domain_envelope(Ê))# 
    end
end
