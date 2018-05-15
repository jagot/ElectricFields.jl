using ElectricFields
using Base.Test

const testfile = joinpath(dirname(@__FILE__), "literate_org_tangled_tests.jl")
if isfile(testfile)
    include(testfile)
else
    error("ElectricFields not properly installed. Please run Pkg.build(\"ElectricFields\") then restart Julia.")
end
