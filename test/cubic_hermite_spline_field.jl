@testset "Cubic Hermite spline fields" begin
    @field(A1d) do
        I₀ = 1.0
        T = 2.0
        σ = 3.0
        σmax = 7.0
    end

    @field(A3d) do
        I₀ = 1.0
        T = 2.0
        σ = 3.0
        σmax = 7.0
        ξ = 0.4
        rotation = (0.3, [1,1,1])
    end

    for A in (A1d, A3d)
        t = timeaxis(A)

        Fv = field_amplitude(A, t)
        Av = vector_potential(A, t)

        B = CubicHermiteSplineField(t, Av, nothing)
        C = CubicHermiteSplineField(t, nothing, Fv)

        # Test that the interpolation works not only at the sample points,
        # but also between them, by making a denser grid.
        tt = range(t[1], stop=t[end], length=2length(t))

        for ttt in (tt, timeaxis(A))
            test_approx_eq(field_amplitude(A, ttt), field_amplitude(B, ttt), rtol=4e-6)
            test_approx_eq(vector_potential(A, ttt), vector_potential(B, ttt), rtol=4e-6)

            test_approx_eq(field_amplitude(A, ttt), field_amplitude(C, ttt), rtol=4e-6)
            test_approx_eq(vector_potential(A, ttt), vector_potential(C, ttt), rtol=4e-6)
        end

        if dimensions(B) == 1
            withenv("UNITFUL_FANCY_EXPONENTS" => false) do
                @test pretty_print_object(B) == """
                                                2200-sample, 1-component Cubic Hermite spline field, t ∈ -22.0 .. 22.0
                                                  - I₀ = 9.9999e-01 = 35.0943 PW cm^-2
                                                    - E₀ = 1.0000e+00 = 514.2192 GV m^-1
                                                    - A₀ = 3.1652e-01
                                                  - ω₀ = 3.1593e+00 Ha = 85.9697 eV (T = 48.1061 as, λ = 14.4219 nm, f = 20.7874 PHz)
                                                  - Intensity FWHM = 7.0832e+00 jiffies = 171.3349 as
                                                  - and a bandwidth of 0.3925 Ha = 10.6799 eV ⟺ 2.5824 PHz ⟺ 33.8563 Bohr = 1.7916 nm ⟹ time–bandwidth product = 4.4245e-01
                                                  - Uₚ = 0.0250 Ha = 679.7651 meV => α = 0.1019 Bohr = 5.3926 pm"""
            end
        else
            withenv("UNITFUL_FANCY_EXPONENTS" => false) do
                @test pretty_print_object(B) == """
                                                2200-sample, 3-component Cubic Hermite spline field, t ∈ -22.0 .. 22.0
                                                  - I₀ = 9.9999e-01 = 35.0943 PW cm^-2
                                                    - E₀ = 1.0000e+00 = 514.2192 GV m^-1
                                                    - A₀ = 3.1652e-01
                                                  - ω₀ = 3.1593e+00 Ha = 85.9697 eV (T = 48.1061 as, λ = 14.4219 nm, f = 20.7874 PHz)
                                                  - Intensity FWHM = 7.0832e+00 jiffies = 171.3349 as
                                                  - and a bandwidth of 0.3925 Ha = 10.6799 eV ⟺ 2.5824 PHz ⟺ 33.8563 Bohr = 1.7916 nm ⟹ time–bandwidth product = 4.4245e-01
                                                  - Uₚ = 0.0069 Ha = 186.7233 meV => α = 0.0890 Bohr = 4.7093 pm"""
            end
        end
    end
end
