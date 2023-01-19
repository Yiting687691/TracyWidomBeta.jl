function step_for(step,A,B,c0,c1,c2,c3,c4,c5)
    TW_cdf=zeros(xl,1)
    TW_pdf=zeros(xl,1)
    TW_cdf[1]=1
    TW_pdf[1]=0
    Δx=xs[2]-xs[1]   
    if step=="trapz"
        for j=1:xl-1
            rhs=(I+(Δx/2)*(A+xs[j]*B))*c0
            lhs=I-(Δx/2)*(A+xs[j+1]*B)
            c0=lhs\rhs
            c0_pdf=(A+xs[j+1]*B)*c0
            TW_cdf[j+1]=c0[end]
            TW_pdf[j+1]=c0_pdf[end]
        end
    elseif step=="bdf3"
        RR = (A+(x0*B))*(6Δx)
        IRR = 11I - RR
        ΔRR = (6Δx^2)*B
        for i = 2:xl
            temp = 18c0 - 9c1 + 2c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            TW_cdf[i] = int(trig_interp(l*pi,c0))(pi)|> real
            TW_pdf[i] = int(trig_interp(l*pi,(RR*c0)/(6Δx)))(pi)|> real
        end
    elseif step=="bdf4"
        RR = (A+(x0*B))*(12Δx)
        IRR = 25I - RR
        ΔRR = (12Δx^2)*B
        for i = 2:xl
            temp = 48c0 - 36c1 + 16c2 - 3c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            TW_cdf[i] = int(trig_interp(l*pi,c0))(pi)|> real
            TW_pdf[i] = int(trig_interp(l*pi,(RR*c0)/(12Δx)))(pi)|> real
        end
    elseif step=="bdf5"
        RR = (A+(x0*B))*(60Δx)
        IRR = 137I - RR
        ΔRR = (60Δx^2)*B
        for i = 2:xl
            temp = 300c0 - 300c1 + 200c2 - 75c3 + 12c4
            c4 = c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            TW_cdf[i] = int(trig_interp(l*pi,c0))(pi)|> real
            TW_pdf[i] = int(trig_interp(l*pi,(RR*c0)/(60Δx)))(pi)|> real
        end
    elseif step=="bdf6"
        RR = (A+(x0*B))*(60Δx)
        IRR = 147I - RR
        ΔRR = (60Δx^2)*B
        for i = 2:xl
            temp = 360c0 - 450c1 + 400c2 - 225c3 + 72c4 - 10c5
            c5 = c4
            c4 = c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            TW_cdf[i] = int(trig_interp(l*pi,c0))(pi)|> real
            TW_pdf[i] = int(trig_interp(l*pi,(RR*c0)/(60Δx)))(pi)|> real
        end
    end
    return TW_cdf,TW_pdf
end
