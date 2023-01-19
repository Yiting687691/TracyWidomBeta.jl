module TracyWidomBeta

using LinearAlgebra
using Statistics
using SparseArrays
using Distributions
using Trapz
using ApproxFun
using SpecialFunctions

struct trig_interp
    L::Float64
    c::Vector{Complex{Float64}}
end

mutable struct int_trig_interp
    L::Float64
    c::Vector{Complex{Float64}}
    C::Complex{Float64}
end

diffvec = (L,m,j) -> ((-floor(m/2):1:floor((m-1)/2))*(1im*pi/L)).^j

function intvec(L,m,j)
    dv = diffvec(L,m,j)
    mm = convert(Int64,floor( m/2 ))
    dv[mm+1] = 1.0
    dv = 1.0./dv
    dv[mm+1] = 0.0
    dv
end

function (tr::trig_interp)(x)
    m = length(tr.c)
    mm = convert(Int64,floor( m/2 ))
    ex = exp.(-1im*pi*mm*x/tr.L)
    ex1 = exp.(1im*pi*x/tr.L)
    sum = tr.c[1]*ex
    for i = 2:length(tr.c)
        ex  =  ex.*ex1
        sum += tr.c[i]*ex
    end
    return sum
end 

function (tr::int_trig_interp)(x)
    trig_interp(tr.L,tr.c)(x) + tr.C*x
end

function int(tr::trig_interp)
    m = length(tr.c)
    mm = convert(Int64,floor( m/2 ))
    tr = int_trig_interp(tr.L,intvec(tr.L,m,1).*tr.c,tr.c[mm+1])
    tr.c[mm+1] = -tr(0.0)
    tr
end

function int(tr::trig_interp,x)
    m = length(tr.c)
    mm = convert(Int64,floor( m/2 ))
    ex = exp.(-1im*pi*mm*x/tr.L)
    ex1 = exp.(1im*pi*x/tr.L)
    sum = tr.c[1]*(ex - 1)*(x/log.(ex))
    for i = 2:length(tr.c)
        ex  =  ex.*ex1
        sum += abs(ex - 1.0) > 1e-14 ? tr.c[i]*(ex - 1)*(x/log.(ex)) : 0.0
    end
    return sum + x*tr.c[mm+1]
end

include("initial_gen.jl")
include("matrix_gen.jl")
include("TW.jl")

export TW

end
