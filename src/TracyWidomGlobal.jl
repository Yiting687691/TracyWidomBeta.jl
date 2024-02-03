function θ_res(β,x0)
    
    Φ = x -> erf.(x/sqrt(2))/2 + 0.5
    g = (t) -> Φ( (x0 - cot(t).^2)./sqrt.(4/β * cot.(t)) )
    h0 = (t) -> t < 0 ? 0.0 : (0 <= t < π/2 ? g(t) : 1.0)
    domain = 0..π
    cheb_space = Chebyshev(domain)
    cheb_piecewise_gau = Fun(t -> h0(t), cheb_space)
    n_θ = length(cheb_piecewise_gau.coefficients) + 20
    cheb_piecewise_gau = Fun(t -> h0(t), cheb_space, n_θ)
    return n_θ, cheb_piecewise_gau.coefficients
    
end

function odo(β,x0,xe,n_x,θ0,θe,n_θ,D_1,S_01,dz)
    
   # C2
    θ = Fun()
    D_2 = (D_1)^2
    M1_θ2 = Multiplication((2/β)*sin((θe-θ0)*(θ+1)+2*θ0)*(sin((θe-θ0)*(θ+1)/2+θ0))^2-(cos((θe-θ0)*(θ+1)/2+θ0))^2, 
        Ultraspherical(1))
    m_C2 = bandwidths(M1_θ2[1:2,1:2])[1]
    M2_θ2 = Multiplication((sin((θe-θ0)*(θ+1)/2+θ0))^4, Ultraspherical(2))
    m_C22 = bandwidths(M2_θ2[1:2,1:2])[1]
    S_12 = Conversion(Ultraspherical(1), Ultraspherical(2))
    C_2 = sparse((8/(β*(θe-θ0)^2))*M2_θ2[1:n_θ,1:n_θ+m_C22]*D_2[1:n_θ+m_C22,1:n_θ]+(2/(θe-θ0))*S_12[1:n_θ,1:n_θ+m_C2]
        *M1_θ2[1:n_θ+m_C2,1:n_θ+m_C2]*D_1[1:n_θ+m_C2,1:n_θ])
    C_2 = convert(SparseMatrixCSC{Float64, Int64},dropzeros(dz.(C_2)))

    # A3
    x = Fun()
    M0_x = Multiplication((x0-xe)*(x+1)/2+xe, Chebyshev())
    m_A3 = bandwidths(M0_x[1:2,1:2])[1]
    A_3 = sparse(S_01[1:n_x,1:n_x+m_A3]*M0_x[1:n_x+m_A3,1:n_x])
    A_3 = convert(SparseMatrixCSC{Float64, Int64},dropzeros(dz.(A_3)))

    # C3
    M1_θ3 = Multiplication((sin((θe-θ0)*(θ+1)/2+θ0))^2, Ultraspherical(1))
    m_C3 = bandwidths(M1_θ3[1:2,1:2])[1]
    C_3 = sparse((2/(θe-θ0))*S_12[1:n_θ,1:n_θ+m_C3]*M1_θ3[1:n_θ+m_C3,1:n_θ+m_C3]*D_1[1:n_θ+m_C3,1:n_θ])
    C_3 = convert(SparseMatrixCSC{Float64, Int64},dropzeros(dz.(C_3)))
    
    return C_2, A_3, C_3
    
end

function odo_up(A_1,C_1,A_2,C_2,A_3,C_3,F,B_x,K,B_θ,G)
    
    C̄_2 = C_2-C_2[:,1]*B_θ
    C̃_2 = C̄_2[1:end-1,2:end]

    Ā_3 = A_3-A_3[:,1]*B_x
    Ã_3 = Ā_3[1:end-1,2:end]

    C̄_3 = C_3-C_3[:,1]*B_θ
    C̃_3 = C̄_3[1:end-1,2:end]

    F̄ = F-(A_1[:,1]*K*C_1'+A_2[:,1]*K*C_2'+A_3[:,1]*K*C_3')-((A_1-A_1[:,1]*B_x)*G'*
        (C_1[:,1])'+(A_2-A_2[:,1]*B_x)*G'*(C_2[:,1])'+(A_3-A_3[:,1]*B_x)*G'*(C_3[:,1])')
    F̃ = F̄[1:end-1,1:end-1]
    
    return C̃_2, Ã_3, C̃_3, F̃

end

function recover(xn, B̃_x, K, B̃_θ, G)
    
    K_2=K[:,2:end]
    K_1=K[:,1]
    G_2=G[:,2:end]
    G_1=G[:,1]
    X_12=K_2-B̃_x*xn
    X_21=G_2'-xn*B̃_θ'
    X_11=G_1'-X_12*B̃_θ'
    top_row = hcat(X_11, X_12)
    bottom_rows = hcat(X_21, xn)
    X=vcat(top_row, bottom_rows)
    return X
    
end
