using StaticArrays

if VERSION < v"1.6.0"
    sincospi(x) = (sinpi(x), cospi(x))
end

@testset "Analytic continuation" begin
    @field(A) do
        I₀ = 2.0
        λ = 100.0
        τ = 3.0
        σmax = 6.0
    end

    @field(A2) do
        I₀ = 2.0
        λ = 100.0
        τ = 3.0
        σmax = 6.0
        rotation = π/4, [0,0,1]
    end

    @field(B) do
        I₀ = 2.0
        λ = 100.0
        τ = 3.0
        toff = 2.5
        tmax = 4.0
        env = :trunc_gauss
    end

    @field(C) do
        I₀ = 2.0
        λ = 100.0
        cycles = 5.0
        env = :cos²
    end

    @field(D) do
        I₀ = 2.0
        λ = 100.0
        ramp = 2.0
        flat = 3.0
        env = :trapezoidal
    end

    for t = (0.4 + 1.0im, 0.4 - 1.0im)
        for f in (A,B)
            α = f.env.α
            ω = photon_energy(A)

            At = f.A₀*exp(-α*t^2)*sin(ω*t)
            @test vector_potential(f, t) ≈ At

            Ft = -f.A₀*exp(-α*t^2)*(-α*2t*sin(ω*t) + ω*cos(ω*t))
            @test field_amplitude(f, t) ≈ Ft
        end

        # Test that complex derivatives works for analytic
        # continuation of vector-valued fields.
        Av = vector_potential(A2, t)
        @test Av isa SVector{3}
        @test Av == [0, 0, vector_potential(A, t)]
        Fv = field_amplitude(A2, t)
        @test Fv isa SVector{3}
        @test Fv == [0, 0, field_amplitude(A, t)]

        let ω = photon_energy(C)
            nT = C.env.cycles*C.env.period
            s,c = sincospi(t/nT)

            At = C.A₀*c^2*sin(ω*t)
            @test vector_potential(C, t) ≈ At

            Ft = -C.A₀*c*(-2π/nT*s*sin(ω*t) + c*ω*cos(ω*t))
            @test field_amplitude(C, t) ≈ Ft
        end

        AC = A + C
        @test vector_potential(AC, t) ≈ vector_potential(A, t) + vector_potential(C, t)
        @test field_amplitude(AC, t) ≈ field_amplitude(A, t) + field_amplitude(C, t)
    end

    let ω = photon_energy(C)
        let t = 2.5+0.4im
            At = D.A₀*sin(ω*t)
            @test vector_potential(D, t) ≈ At

            Ft = -D.A₀*ω*cos(ω*t)
            @test field_amplitude(D, t) ≈ Ft
        end
    end
end
