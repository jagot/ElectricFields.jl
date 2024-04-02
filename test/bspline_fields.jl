@testset "B-spline fields" begin
    @testset "$(n)d fields" for n in (1,3)
        @field(F) do
            ω = 1.0
            I₀ = 1.0
            σ = 7.0
            σoff = 4.0
            σmax = 6.0
            env = :trunc_gauss
        end
        if n == 3
            F = rotate(F, ElectricFields.compute_rotation((π/3, [0.4,1,0])))
        end

        t = timeaxis(F)
        Av = vector_potential(F, t)
        Fv = field_amplitude(F, t)

        spline_order = 7
        num_steps = length(t)
        num_knots = 100

        BB = ElectricFields.BSpline(ElectricFields.LinearKnotSet(spline_order, t[1], t[end], num_knots))
        @testset "B-spline $(case)" for case in ("full range", "restricted")
            B = if case == "restricted"
                BB[:,begin+1:end-1]
            else
                BB
            end

            FB = BSplineField(B, t, Av)

            @test dimensions(FB) == n

            @test pretty_print_object(FB) ==
                  if case == "restricted"
                  "B-spline field expanded over\n  Restriction to basis functions 2:105 of BSpline basis with LinearKnotSet(Float64) of order k = 7 on -42.0 .. 42.0 (100 intervals)"
                  else
                      "B-spline field expanded over\n  BSpline basis with LinearKnotSet(Float64) of order k = 7 on -42.0 .. 42.0 (100 intervals)"
                  end

            Arec = vector_potential(FB, t)
            Frec = field_amplitude(FB, t)

            test_approx_eq(Av, Arec, rtol=1e-4)
            test_approx_eq(Fv, Frec, rtol=1e-4)
        end
    end
end
