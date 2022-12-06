function TW(β;method="finite",type="cdf")
    if method=="finite" && type=="cdf"
        sol=finite_dis(β,1000)
        TW=sol[1]
    elseif method=="finite" && type=="pdf"
        sol=finite_dis(β,1000)
        TW=sol[2]
    elseif method=="spectral" && type=="cdf"
        sol=spectral_dis(β,1000)
        TW=sol[1]
    elseif method=="spectral" && type=="pdf"
        sol=spectral_dis(β,1000)
        TW=sol[2]
    else
        TW="Input valid method"
    end
    return TW
end
