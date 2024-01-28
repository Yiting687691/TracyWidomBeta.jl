module TracyWidomBeta

using LinearAlgebra, Statistics, SparseArrays, Distributions, Trapz, ApproxFun, SpecialFunctions, FFTW

include("time_gen.jl")
include("initial_gen.jl")
include("matrix_gen.jl")
include("step_for.jl")
include("Fourier_interp.jl")
include("TW.jl")

export TW

end
