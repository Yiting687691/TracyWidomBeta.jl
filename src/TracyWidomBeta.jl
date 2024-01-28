module TracyWidomBeta

using LinearAlgebra, Statistics, SparseArrays, Distributions, Trapz, ApproxFun, SpecialFunctions, FFTW

include("initial_gen.jl")
include("matrix_gen.jl")
include("step_for.jl")
include("Fourier_interp.jl")
include("TW.jl")

export TW

function time_gen(x1,x2,Δx)
    xs=x1:Δx:x2 |> Array
    if abs(xs[end])<abs(x2)
        xs=x1:Δx:(x2+Δx) |> Array
    end
    return xs
end

end
