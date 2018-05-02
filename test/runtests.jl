using ElectricFields
using Base.Test

const testfile = joinpath(dirname(@__FILE__), "..", "deps", "build", "tests.jl")
if isfile(testfile)
    include(testfile)
else
    error("ElectricFields not properly installed. Please run Pkg.build(\"ElectricFields\") then restart Julia.")
end
