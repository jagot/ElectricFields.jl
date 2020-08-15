var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ElectricFields","category":"page"},{"location":"#ElectricFields","page":"Home","title":"ElectricFields","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The idea of this package is to provide an interface between the reality and calculations. In the calculations, it is useful to represent fields in terms of cycles of a fundamental frequency, which yields a timebase. E.g. one might use laser pulses of 800 nm, in an experiment, which has a period time of about 2.66 fs. It is, however, easier to calculate in normalized time, and relate all other quantities of interest (such as ionization potential, &c) to this time scale.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package provides a simple DSL that requires just a handful of parameters, that can be given in any unit system (thanks to Unitful.jl). Different fields can be combined in any way that is physically reasonable, to recreate complicated experimental situations. Everything can then be converted to normalized time, for use in calculations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ElectricFields]","category":"page"},{"location":"#ElectricFields.BK7","page":"Home","title":"ElectricFields.BK7","text":"BK7\n\nBorosilicate glass is commonly used in optical lenses\n\n\n\n\n\n","category":"constant"},{"location":"#ElectricFields.SiO₂","page":"Home","title":"ElectricFields.SiO₂","text":"SiO₂\n\nFused silica\n\n\n\n\n\n","category":"constant"},{"location":"#ElectricFields.GaussianEnvelope","page":"Home","title":"ElectricFields.GaussianEnvelope","text":"GaussianEnvelope\n\nA Gaussian pulse is given by\n\nI_0expleft(-fract^22sigma^2right)\n\nwhere the standard deviation σ is related to the FWHM duration τ of the intensity envelope as\n\nsigma = fractau2sqrt2ln 2\n\nFurthermore, the amplitude standard deviation σ is proportional to the intensity ditto: sigma = sqrt2sigma. Therefore, the amplitude envelope is given by\n\nE_0expleft(-fract^22sigma^2right)\n=E_0expleft(-fract^24sigma^2right)\n=E_0expleft(-frac2ln2t^2tau^2right)\n\nSince a Gaussian never ends, we specify how many σ we require; the resulting time window will be rounded up to an integer amount of cycles of the fundamental.\n\n\n\n\n\n","category":"type"},{"location":"#ElectricFields.Medium","page":"Home","title":"ElectricFields.Medium","text":"Medium(B, C)\n\nThe Sellmeier equations are used to describe dispersion in glass using a series of resonances:\n\nn^2(lambda) =\n1 + sum_i frac B_ilambda^2lambda^2-C_i\n\n\n\n\n\n","category":"type"},{"location":"#ElectricFields.PaddedField","page":"Home","title":"ElectricFields.PaddedField","text":"  PaddedField(field, a, b)\n\nWrapper around any electric field, padded with a units of time before, and b units of time after the ordinary span of the field.\n\n\n\n\n\n","category":"type"},{"location":"#ElectricFields.WrappedField","page":"Home","title":"ElectricFields.WrappedField","text":"WrappedField(field, a, b)\n\nWrapper around any electric field\n\n\n\n\n\n","category":"type"},{"location":"#ElectricFields.calc_params!-Tuple{Dict{Symbol,Any}}","page":"Home","title":"ElectricFields.calc_params!","text":"calc_params!(field_params)\n\nThis function performs the calculation of different quantities from the information provided.\n\nThe ponderomotive potential U_p is the cycle-average quiver energy of a free electron in an electromagnetic field. It is given by\n\nU_p =\nfrace^2E_0^24momega^2=frac2e^2cvarepsilon_0mtimesfracI4omega^2\n\nor, in atomic units,\n\nU_p = fracI4omega^2\n\n\n\n\n\n","category":"method"},{"location":"#ElectricFields.dispersion-Union{Tuple{F}, Tuple{Medium,Union{Unitful.Quantity{T,𝐋,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐋,U}} where S where L} where U where T,AbstractArray{F,1}}, Tuple{Medium,Union{Unitful.Quantity{T,𝐋,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐋,U}} where S where L} where U where T,AbstractArray{F,1},Union{Unitful.Quantity{T,𝐓^-1,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐓^-1,U}} where S where L} where U where T}} where F<:(Union{Unitful.Quantity{T,𝐓^-1,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐓^-1,U}} where S where L} where U where T)","page":"Home","title":"ElectricFields.dispersion","text":"dispersion(m, d, f[, f₀=0u\"Hz])\n\nCalculate dispersion through a medium m of length d. Optionally, remove central frequency k-vector, to keep pulse temporally centred.\n\n\n\n\n\n","category":"method"},{"location":"#ElectricFields.keldysh-Tuple{ElectricFields.AbstractField,Union{Unitful.Quantity{T,𝐋^2 𝐌 𝐓^-2,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐋^2 𝐌 𝐓^-2,U}} where S where L} where U where T}","page":"Home","title":"ElectricFields.keldysh","text":"keldysh(f, Iₚ)\n\nThe Keldysh parameter relates the strength of a dynamic electric field to that of the binding potential of an atom. It is given by\n\ngamma = sqrtfracI_p2U_p\n\nwhere I_p is the ionization potential of the atom and U_p is the ponderomotive potential of the dynamic field.\n\n\n\n\n\n","category":"method"},{"location":"#ElectricFields.spectrum-Tuple{ElectricFields.GaussianEnvelope}","page":"Home","title":"ElectricFields.spectrum","text":"spectrum(env::GaussianEnvelope)\n\nGaussians belong to the Schwartz class, i.e. functions who, under Fourier transform, are mapped back to the same space. That is to say, the Fourier transform of a Gaussian is a Gaussian:\n\nexp(-alpha t^2) leftrightarrow\nfrac1sqrt2alpha\nexpleft(-fracomega^24alpharight)\n\nComparing with the above, we find that the spectral standard deviation\n\nOmega = sqrt2alpha = frac2sqrtln 2tau\n\nand the Gaussian function in the spectral domain is thus\n\nE(omega) =\nfracE_0tau2sqrtln 2\nexpleft-frac(omegatau)^28ln2right\n\n\n\n\n\n","category":"method"},{"location":"#ElectricFields.test_field_parameters-Tuple{Any,Any}","page":"Home","title":"ElectricFields.test_field_parameters","text":"test_field_parameters(field_params, set)\n\nThis function ensures that one and only one of \"competing\" quantities is specified.\n\n\n\n\n\n","category":"method"},{"location":"#ElectricFields.@namespace!-Tuple{Any,Any}","page":"Home","title":"ElectricFields.@namespace!","text":"@namespace!(exprs, params)\n\nThis macro uses the dictionary params as a \"namespace\", i.e. all symbols are assumed to be keys in this dictionary. We need to escape the generated expression tree to modify in the scope of the caller, not the global scope.\n\n\n\n\n\n","category":"macro"}]
}