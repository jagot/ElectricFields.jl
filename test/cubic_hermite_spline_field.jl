@testset "Cubic Hermite spline fields" begin
    @field(A) do
        I₀ = 1.0
        T = 2.0
        σ = 3.0
        σmax = 7.0
    end
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

    @test pretty_print_object(A) == """
                                    Linearly polarized field with
                                      - I₀ = 1.0000e+00 au = 3.5094452e16 W cm^-2 =>
                                        - E₀ = 1.0000e+00 au = 514.2207 GV m^-1
                                        - A₀ = 0.3183 au
                                      – a Fixed carrier @ λ = 14.5033 nm (T = 48.3777 as, ω = 3.1416 Ha = 85.4871 eV, f = 20.6707 PHz)
                                      – and a Gaussian envelope of duration 170.8811 as (intensity FWHM; ±7.33σ)
                                      – and a bandwidth of 0.3925 Ha = 10.6797 eV ⟺ 2.5823 PHz ⟺ 34.2390 Bohr = 1.8119 nm
                                      – Uₚ = 0.0253 Ha = 689.2724 meV => α = 0.1013 Bohr = 5.3617 pm"""
end
