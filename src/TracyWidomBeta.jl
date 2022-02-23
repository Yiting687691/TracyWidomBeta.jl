module TracyWidomBeta
using LinearAlgebra,Statistics,SparseArrays,Random,Distributions,Plots,Printf,Trapz,RandomMatrices
export CDF
export PDF
include("CDF.jl")
include("PDF.jl")
end
