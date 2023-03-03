function TW(β;cheb=10^3,method="finite",interp=true,step="trapz",pdf=false,x0=13.0/(sqrt(β)),xN=-10.0,Δx_f=-0.001,Δx_s=-0.001,M_f=Int(floor(-1/(Δx_f))),M_s=8000,l=10)
    
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
