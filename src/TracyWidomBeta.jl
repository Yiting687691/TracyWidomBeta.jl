module TracyWidomBeta

using LinearAlgebra, Statistics, SparseArrays, Distributions, Trapz, ApproxFun, SpecialFunctions, FFTW

include("AuxiliaryFunctions.jl")
include("TracyWidomGlobal.jl")

export TW,TW_glob

function TW(β;cheb=10^3,method="finite",interp=true,step="trapz",pdf=false,x0=floor(13.0/(sqrt(β))),xN=-10.0,Δx_f=-0.001,Δx_s=-0.001,M_f=Int(floor(-1/(Δx_f))),M_s=8000,l=10)
    
    # Set up the time domain
    if method=="finite"
        xs=time_gen(x0,xN,Δx_f)
    elseif method=="spectral"
        xs=time_gen(x0,xN,Δx_s)
    end
    xl=length(xs)
    
    # Set up the initial conditions
    (c0,c1,c2,c3,c4,c5,h,θ)=initial_gen(β;method,x0,Δx_f,Δx_s,M_f,M_s,l)

    # Set up the matrices
    (A,B)=matrix_gen(β;method,M_f,M_s,h,θ,l)

    # Step forward in time
    (TW_cdf,TW_pdf)=step_for(method,step,A,B,x0,M_s,xs,xl,c0,c1,c2,c3,c4,c5,l)

    # Interpolation
    (cdf_cheb,pdf_cheb)=Fourier_interp(xs,cheb,TW_cdf,TW_pdf)
    
    if interp && pdf==false
        return cdf_cheb
    elseif interp && pdf
        return pdf_cheb
    elseif interp==false && pdf==false
        return xs,TW_cdf
    elseif interp==false && pdf
        return xs,TW_pdf
    else
        return "Input valid arguments"
    end
end

function TW_glob(β;x0=round.(((7*β*log(β)+47)^(2/3))*β^(-2/3); digits = 2),xe = -10,θ0 = 0, n_x=Int(ceil(4*β/3+50/3)), B=10, E=1)
    
    dz = (t) -> abs(t) <= 1e-16 ? 0 : t
    cheb_coeffr = zeros(n_x,B)
    cheb_coeffl = zeros(n_x,B)
    n_θ, coeffb = θ_res(β,x0)
    s = (x0 - xe)/B
    endx = [x0 - i*s for i in 0:B]
    endt = [θ0 + i*π for i in 0:E]
    
    # A1
    D_1 = Derivative()
    D_1 = D_1:Chebyshev()
    A_1 = sparse(D_1[1:n_x,1:n_x])
    A_1 = convert(SparseMatrixCSC{Float64, Int64},dropzeros(dz.(A_1)))
    
    # C1
    S_02 = Conversion(Chebyshev(), Ultraspherical(2))
    C_1 = sparse((2/s)*S_02[1:n_θ,1:n_θ])
    C_1 = convert(SparseMatrixCSC{Float64, Int64},dropzeros(dz.(C_1)))
    
    # A2
    S_01 = Conversion(Chebyshev(), Ultraspherical(1))
    A_2 = sparse(S_01[1:n_x,1:n_x])
    A_2 = convert(SparseMatrixCSC{Float64, Int64},dropzeros(dz.(A_2)))
    
    #F
    F = convert(SparseMatrixCSC{Float64, Int64},dropzeros(dz.(sparse(zeros(n_x,n_θ)))))
    
    # Bx
    B_x = ones(1,n_x)
    B_x = convert(SparseMatrixCSC{Float64, Int64},B_x)
    B̃_x = B_x[:,2:end]
    
    # Bθ
    B_θ = ones(1,n_θ)
    B_θ[1,:] = [(-1)^(i-1) for i in 1:n_θ]
    B_θ = convert(SparseMatrixCSC{Float64, Int64},B_θ)
    B̃_θ = B_θ[:,2:end]
    
    # Initial K and G
    K = zeros(1,n_θ)
    G = zeros(1,n_x)
    
    # Ã1
    Ā_1 = A_1-A_1[:,1]*B_x
    Ã_1 = Ā_1[1:end-1,2:end]
    
    # C̃1
    C̄_1 = C_1-C_1[:,1]*B_θ
    C̃_1 = C̄_1[1:end-1,2:end]
    
    # Ã2
    Ā_2 = A_2-A_2[:,1]*B_x
    Ã_2 = Ā_2[1:end-1,2:end]
    
    v = [(-1)^(i-1) for i in 1:n_x]
    
    for i = 1:E
        θ0̂ = endt[i]
        θê = endt[i+1]
        for j = 1:B
            x0̂ = endx[j]
            xê = endx[j+1]
            coefft = coeffb
            coeffl = cheb_coeffl[:,j]
            C_2, A_3, C_3 = odo(β,x0̂,xê,n_x,θ0̂,θê,n_θ,D_1,S_01,dz)
            K[1,:] = coefft
            K = convert(SparseMatrixCSC{Float64, Int64},K)
            G[1,:] = coeffl
            G = convert(SparseMatrixCSC{Float64, Int64},G)
            C̃_2,Ã_3,C̃_3,F̃ = odo_up(A_1,C_1,A_2,C_2,A_3,C_3,F,B_x,K,B_θ,G)
            Kr = kron(C̃_1,Ã_1)+kron(C̃_2,Ã_2)+kron(C̃_3,Ã_3)
            Xb = Matrix(Kr)\vec(Matrix(F̃))
            xn = reshape(Xb,n_x-1,n_θ-1)
            X = recover(xn, B̃_x, K, B̃_θ, G)
            coeffr = vec(sum(X, dims=2))
            cheb_coeffr[:,j] = coeffr
            coeffb = vec(X'*v)
        end
        cheb_coeffl = cheb_coeffr
        coeffb = vec(zeros(n_θ,1))
        coeffb[1] = 1
    end
    TW_cdf = 0
    TW_cdf = sum([Fun(x0-i*s..x0-(i-1)*s, cheb_coeffr[:,i]) for i in 1:B])
    return TW_cdf
    
end

end
