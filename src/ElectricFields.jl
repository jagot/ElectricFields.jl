__precompile__()

module ElectricFields

const codefile = joinpath(dirname(@__FILE__), "..", "deps", "build", "code.jl")
if isfile(codefile)
    include(codefile)
else
    error("ElectricFields not properly installed. Please run Pkg.build(\"ElectricFields\") then restart Julia.")
end

end # module
