module TracyWidomBeta

using LinearAlgebra
using Statistics
using SparseArrays
using Distributions
using Trapz
using ApproxFun
using SpecialFunctions
import Plots

include("one_step.jl")
include("finite_cdf.jl")
include("BDF4_cdf.jl")
include("TW.jl")
include("model_test.jl")
include("density_theta.jl")
include("water_den.jl")
include("water_TW.jl")
include("den_other.jl")

export TW

end
