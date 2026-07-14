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

    test_approx_eq(field_amplitude(A, tt), field_amplitude(B, tt), rtol=4e-6)
    test_approx_eq(vector_potential(A, tt), vector_potential(B, tt), rtol=4e-6)

    test_approx_eq(field_amplitude(A, tt), field_amplitude(C, tt), rtol=4e-6)
    test_approx_eq(vector_potential(A, tt), vector_potential(C, tt), rtol=4e-6)
end
