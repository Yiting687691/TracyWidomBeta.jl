function step_for(method,step,A,B,x0,M_s,xs,xl,c0,c1,c2,c3,c4,c5,l)
    TW_cdf=zeros(xl,1)
    TW_pdf=zeros(xl,1)
    TW_cdf[1]=1
    TW_pdf[1]=0
    Δx=xs[2]-xs[1]
    N=-floor(M_s/2):1:floor((M_s-1)/2)
    integ=zeros(ComplexF64,length(N))
    for j=1:length(N)
        if N[j]==0
            integ[j]=pi
        else
            integ[j]=l/(im*N[j])*(exp(im*N[j]*pi/l)-1)
        end
    end
    if step=="trapz"
        for i=1:xl-1
            rhs=(I+(Δx/2)*(A+xs[i]*B))*c0
            lhs=I-(Δx/2)*(A+xs[i+1]*B)
            c0=lhs\rhs
            c0_pdf=(A+xs[i+1]*B)*c0
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1] = c0[end]
                TW_pdf[i+1] = c0_pdf[end]
            end
        end
    elseif step=="bdf3"
        RR = (A+(x0*B))*(6Δx)
        IRR = 11I - RR
        ΔRR = (6Δx^2)*B
        for i = 1:xl-1
            temp = 18c0 - 9c1 + 2c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            c0_pdf=(A+xs[i+1]*B)*c0
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1]=c0[end]
                TW_pdf[i+1]=c0_pdf[end]
            end
        end
    elseif step=="bdf4" && method=="spectral"
        RR = (A+(x0*B))*(12Δx)
        IRR = 25I - RR
        ΔRR = (12Δx^2)*B
        for i = 1:xl-1
            temp = 48c0 - 36c1 + 16c2 - 3c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            c0_pdf=(A+xs[i+1]*B)*c0
            TW_cdf[i+1] = sum(c0.*integ)|> real
            TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
        end
    elseif step=="bdf5" && method=="spectral"
        RR = (A+(x0*B))*(60Δx)
        IRR = 137I - RR
        ΔRR = (60Δx^2)*B
        for i = 1:xl-1
            temp = 300c0 - 300c1 + 200c2 - 75c3 + 12c4
            c4 = c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            c0_pdf=(A+xs[i+1]*B)*c0
            TW_cdf[i+1] = sum(c0.*integ)|> real
            TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
        end
    elseif step=="bdf6" && method=="spectral"
        RR = (A+(x0*B))*(60Δx)
        IRR = 147I - RR
        ΔRR = (60Δx^2)*B
        for i = 1:xl-1
            temp = 360c0 - 450c1 + 400c2 - 225c3 + 72c4 - 10c5
            c5 = c4
            c4 = c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            c0_pdf=(A+xs[i+1]*B)*c0
            TW_cdf[i+1] = sum(c0.*integ)|> real
            TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
        end
    elseif step=="bdf4" || step=="bdf5" || step=="bdf6" && method=="finite"
        RR = (A+(x0*B))*(6Δx)
        IRR = 11I - RR
        ΔRR = (6Δx^2)*B
        for i = 1:xl-1
            temp = 18c0 - 9c1 + 2c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            c0_pdf=(A+xs[i+1]*B)*c0
            TW_cdf[i+1]=c0[end]
            TW_pdf[i+1]=c0_pdf[end]
        end
        z = "To ensure convergence, BDF3 is used instead."
        println(z)
    end
    return TW_cdf,TW_pdf
end
