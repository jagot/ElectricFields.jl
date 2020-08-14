@testset "Parametric field creation" begin
    #    Test parametric creation of fields, both from inside a function and
    #    from global scope.

    function test_parametric_fields(l)
        Is = range(0u"W/cm^2",stop=1e6u"W/cm^2", length=l)

        for (i,I₀) in enumerate(Is)
            @field(IR) do
                λ = 800.0u"nm"
                I₀ = I₀
                τ = 6.2u"fs"
                Tmax = 5
            end
            @test intensity(IR) == I₀
        end
    end
    test_parametric_fields(11)
    test_parametric_fields(8)

    Is = range(0u"W/cm^2",stop=1e6u"W/cm^2", length=11)

    for (i,I₀) in enumerate(Is)
        @field(IR) do
            λ = 800.0u"nm"
            I₀ = I₀
            τ = 6.2u"fs"
            Tmax = 5
        end
        @test intensity(IR) == I₀
    end

    @field(IR) do
        λ = 800.0u"nm"
        I₀ = 1e14u"W/cm^2"
        τ = 6.2u"fs"
        Tmax = 5
    end

    @test dimensions(IR) == 2
end
