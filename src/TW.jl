function TW(beta;method="finite")
    if method=="finite"
        sol=finite_cdf(beta)
        TW=sol[1]
    elseif method=="BDF4"
        sol=BDF4_cdf(beta)
        TW=sol[1]
    else
        TW="Input valid method"
    end
    return TW
end
