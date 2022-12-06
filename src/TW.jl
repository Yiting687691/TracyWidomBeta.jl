function TW(β,cheb_coef;method="finite")
    if method=="finite"
        sol=finite_dis(β,cheb_coef)
        TW=sol
    elseif method=="spectral"
        sol=spectral_dis(β,cheb_coef)
        TW=sol
    else
        TW="Input valid method"
    end
    return TW
end
