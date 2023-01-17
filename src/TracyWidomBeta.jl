module TracyWidomBeta

using LinearAlgebra
using Statistics
using SparseArrays
using Distributions
using Trapz
using ApproxFun
using SpecialFunctions

include("initial_gen.jl")
include("matrix_gen.jl")
include("TW.jl")

export TW

end
