import ElectricFields: oneaunit,
    get_unitful_quantity,
    shift_unit, si_round,
    au2si, au2si_round,
    I2si_round, Iaustrip

@testset "Units" begin
    @test oneaunit(800u"nm") == u"bohr"
    @test oneaunit(1e14u"W/cm^2") == ElectricFields.Iau

    params = Dict{Symbol,Any}(:τ => 6.2u"fs",:λ => 1.0)
    @test get_unitful_quantity(params, :τ) == 6.2u"fs"
    @test get_unitful_quantity(params, :λ) == 1.0u"bohr"

    @test shift_unit(u"nm", -2) == (u"pm", -1)
    @test shift_unit(u"nm", -3) == (u"pm", -1)
    @test shift_unit(u"nm", -4) == (u"fm", -2)
    @test shift_unit(u"nm", 3) == (u"μm", 1)
    @test shift_unit(u"nm", 1) == (u"nm", 0)
    @test shift_unit(u"ym", -30) == (u"ym", 0)

    @test au2si(1, u"nm") ≈ 0.0529177210903u"nm"
    @test au2si(3, u"bohr") == 3u"bohr"
    @test au2si(u"bohr") == 1u"bohr"
    @test au2si_round(1, u"nm") == "52.9177 pm"

    @test I2si_round(1.0) == "35.0945 PW cm⁻²"
    @test I2si_round(Iaustrip(1e14u"W/cm^2")) == "100.0000 TW cm⁻²"
    
    @test Iaustrip(ElectricFields.Iau) ≈ 1.0 rtol=1e-14
end
