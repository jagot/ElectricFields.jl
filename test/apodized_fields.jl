@testset "Apodized fields" begin
    @field(F) do
        ω = 1.0
        I₀ = 1.0
        ramp = 0.0
        flat = 3.0
        env = :tophat
    end

    t = timeaxis(F)
    tplot = ustrip.(auconvert.(u"fs", t))

    tmin = auconvert(u"fs", 1.0)
    tmax = auconvert(u"fs", 14.0)

    x = ElectricFields.winx.(austrip(tmin), austrip(tmax), austrip.(t))

    pre_sel = findall(<(1.0), t)
    mid_sel = findall(∈(1.0 .. 14.0), t)
    mid_sel1 = findall(∈(1.0 .. 7.5), t)
    mid_sel2 = findall(∈(7.5 .. 14.0), t)
    post_sel = findall(>(14.0), t)

    ElectricFields.@cosine_sum_window Zero () "Zero"
    ElectricFields.@cosine_sum_window One (1.0,) "One"

    @testset "Window: $(window)" for window in (ElectricFields.Hann(), ElectricFields.Hamming(),
                                                ElectricFields.Blackman(),
                                                ElectricFields.BlackmanExact(),
                                                ElectricFields.Nuttall(),
                                                ElectricFields.BlackmanNuttall(),
                                                ElectricFields.BlackmanHarris(),
                                                ElectricFields.Kaiser(3), ElectricFields.Kaiser(2),
                                                ElectricFields.Rect(), Zero(), One())
        w = ElectricFields.window_value.(window, x)
        ∂w = ElectricFields.window_derivative.(window, x)

        @test all(iszero, w[pre_sel])
        @test all(iszero, ∂w[pre_sel])

        if window isa Union{ElectricFields.Rect,One}
            @test all(isone, w[mid_sel])
            @test all(iszero, ∂w[mid_sel])
        elseif window isa Zero
            @test all(iszero, w[mid_sel])
            @test all(iszero, ∂w[mid_sel])
        else
            @test all(<(1), w[mid_sel])
            # This is of course not a proper test of the accuracy of
            # the derivative, just that the behaviour is reasonable.
            @test all(≥(0), ∂w[mid_sel1])
            @test all(≤(0), ∂w[mid_sel2])
        end

        @test all(iszero, w[post_sel])
        @test all(iszero, ∂w[post_sel])
    end

    Fw = ApodizedField(F, tmin, tmax, ElectricFields.Kaiser(3))

    tw = timeaxis(Fw)

    @test parent(Fw) == F
    @test span(Fw) ≈ 1.0 .. 14.0
    @test tw ⊆ t

    Fv = field_amplitude(F, t)
    Av = vector_potential(F, t)

    Fwv = field_amplitude(Fw, t)
    Awv = vector_potential(Fw, t)

    Iv = intensity(F, t)
    Iwv = intensity(Fw, t)

    w = ElectricFields.window_value.(Fw.window, 1.0, 14.0, t)

    @test all(iszero, Fwv[pre_sel])
    # There is no guarantee that the field amplitude is always smaller
    # for the apodized field, due to the derivative term.
    @test all(abs.(Iwv[mid_sel]) .≤ abs.(Iv[mid_sel]))
    # For the same reason, the intensity is not exactly the original
    # intensity (flat for this particular field) times the window
    # value squared, hence the high tolerance.
    @test abs2.(w[mid_sel]) ≈ Iwv[mid_sel] rtol=5e-2
    @test all(iszero, Fwv[post_sel])

    @test all(iszero, Awv[pre_sel])
    @test all(abs.(Awv[mid_sel]) .≤ abs.(Av[mid_sel]))
    @test all(iszero, Awv[post_sel])

    withenv("UNITFUL_FANCY_EXPONENTS" => true) do
        @test string(Fw) == """
        Kaiser(α = 3) window from 1.0000 jiffies = 24.1888 as to 14.0000 jiffies = 338.6438 as of
        Linearly polarized field with
          - I₀ = 1.0000e+00 au = 3.5094452e16 W cm⁻² =>
            - E₀ = 1.0000e+00 au = 514.2207 GV m⁻¹
            - A₀ = 1.0000 au
          – a Fixed carrier @ λ = 45.5634 nm (T = 151.9830 as, ω = 1.0000 Ha = 27.2114 eV, f = 6.5797 PHz)
          – and a /0‾3‾0\\ cycles trapezoidal envelope
          – and a bandwidth of Inf Ha = Inf eV ⟺ Inf Hz ⟺ Inf Bohr = Inf m
          – Uₚ = 0.2500 Ha = 6.8028 eV => α = 1.0000 Bohr = 52.9177 pm"""
    end

    R = ElectricFields.compute_rotation((π/2, [0,1,0]))
    RFw = rotate(Fw, R)
    RFwv = field_amplitude(RFw, t)
    @test RFwv ≈ hcat(Fwv, zeros(length(t), 2))
end
