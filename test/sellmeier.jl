@testset "Sellmeier" begin
    @testset "String representation" begin
        Bstr = string(BK7.B)
        Cstr = join(string.(BK7.C), ", ")
        @test string(BK7) == "Sellmeier(0.0, $(Bstr), [2, 2, 2], [$(Cstr)], [], $(Int)[])"
    end

    @test BK7(0.5876u"μm") ≈ 1.5168 atol=1e-3
    @test SiO₂(0.5876u"μm") ≈ 1.4585 atol=1e-3

    Cn² = ElectricFields.n².(Ref(Calcite), [1064,800,632.8,532,400,355,266]*u"nm")
    Cno = .√(first.(Cn²))
    Cne = .√(last.(Cn²))
    test_approx_eq(Cno, [1.6423,1.6487,1.6557,1.6629,1.6823,1.6951,1.7497], rtol=5e-5)
    test_approx_eq(Cne, [1.4797,1.4821,1.4852,1.4886,1.4974,1.5032,1.5259], rtol=5e-5)

    Qn² = ElectricFields.n².(Ref(Quartz), [1064,800,632.8,532,400,355,266]*u"nm")
    Qno = .√(first.(Qn²))
    Qne = .√(last.(Qn²))
    test_approx_eq(Qno, [1.5341,1.5384,1.5427,1.5469,1.5577,1.5646,1.5916], rtol=5e-5)
    test_approx_eq(Qne, [1.5428,1.5473,1.5517,1.5561,1.5673,1.5744,1.6024], rtol=5e-5)

    Kn² = ElectricFields.n².(Ref(KTP), [1064,532]*u"nm")
    Knx = real(.√(ElectricFields.maybe_complex.(first.(Kn²))))
    Kny = real(.√(ElectricFields.maybe_complex.((e -> e[2]).(Kn²))))
    Knz = real(.√(ElectricFields.maybe_complex.(last.(Kn²))))
    test_approx_eq(Knx, [1.7377,1.7780], rtol=5e-3)
    test_approx_eq(Kny, [1.7453,1.7886], rtol=5e-3)
    test_approx_eq(Knz, [1.8297,1.8887], rtol=5e-3)
end
