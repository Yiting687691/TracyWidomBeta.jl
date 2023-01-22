module TracyWidomBeta

using LinearAlgebra
using Statistics
using SparseArrays
using Distributions
using Trapz
using ApproxFun
using SpecialFunctions

include("time_gen.jl")
include("initial_gen.jl")
include("matrix_gen.jl")
include("step_for.jl")
include("Fourier_interp.jl")
include("TW.jl")

export TW

end
