function Fourier_interp(TW_cdf,TW_pdf)  
    ϕ = x -> (erf.(x) .+ 1.0)/2
    j=PeriodicSegment(xs[end],xs[1])
    S=Laurent(j)
    final_pdf=TW_pdf[2:end] |> reverse
    final_cdf=TW_cdf[2:end] |> reverse
    xx=xs[2:end] |> reverse
    f_pdf=Fun(S,ApproxFun.transform(S,final_pdf))
    f_cdf=Fun(S,ApproxFun.transform(S,final_cdf-ϕ(xx)))
    S2=xs[end]..xs[1] |> Chebyshev
    pdf_cheb=Fun(f_pdf,S2,cheb) |> real
    cdf_cheb=Fun(f_cdf,S2,cheb) + Fun(ϕ,S2) |> real
    return cdf_cheb,pdf_cheb
end
