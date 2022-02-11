module TracyWidomBeta
using LinearAlgebra,Statistics,SparseArrays,Random,Distributions,Plots,Printf,Trapz,RandomMatrices
export F
println("What is the value of beta?")
num1=readline()
beta=num = parse(Float64, num1)
println("What is the value of s?")
num2=readline()
s=parse(Float64, num2)
println("Do you want high accuracy or low accuracy?")
answer=readline()
if answer=="low"
    include("finite.jl")
elseif answer=="high"
    include("BDF4.jl")
else
    println("Try to enter either low or high.")
end

end
