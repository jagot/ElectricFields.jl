import ElectricFields: @namespace!

@testset "Namespace macro" begin
    params = Dict{Symbol,Rational}()

    @namespace!(params) do
        x = 5
        y = 3x
        z = x//3
        if x
            w = x^2
        end
    end

    @test params[:x] == 5
    @test params[:y] == 15
    @test params[:z] == 5//3
    @test params[:w] == 25

    params2 = Dict{Symbol,Quantity}(:λ => 1u"m")

    @namespace!(params2) do
        if λ
            ν = 1/λ
        end
    end

    @test params2[:ν] == 1u"m"^-1

    params3 = Dict{Symbol,Rational}()

    factor = 4
    @namespace!(params3) do
        x = 5
        y = factor*x
    end

    @test params3[:y] == 20

    for i in 1:8
        params4 = Dict{Symbol,Rational}()
        @namespace!(params4) do
            x = 5
            y = i*x
        end

        @test params4[:y] == 5i
    end
end
