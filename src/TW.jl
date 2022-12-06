function TW(β,cheb_coef;method="finite",x0=13.0,xN=-10.0,Δx_f=-0.01,Δx_s=-0.002,M_f=10^3,M_s=8000,l=20)
    if method=="finite"
        TW=finite_dis(β,cheb_coef;x0,xN,Δx_f,M_f)
    elseif method=="spectral"
        TW=spectral_dis(β,cheb_coef;x0,xN,Δx_s,M_s,l)
    else
        TW="Input valid method"
    end
    return TW
end
