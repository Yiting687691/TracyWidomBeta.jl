function step_for(method,step,A,B,x0,M_s,xs,xl,c0,c1,c2,c3,c4,c5,l;x33=0.4358665215,α=(2/sqrt(3))*cos(π/18),x34=1.06858)
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
    elseif step=="TRBDF2"    #TR-BDF2
        for i=1:xl-1
            rhs=(I+(Δx/4)*(A+xs[i]*B))*c0
            lhs=I-(Δx/4)*(A+(xs[i]+(Δx/2))*B)
            temp=lhs\rhs
            rhs=(4/3)*temp-(1/3)*c0
            lhs=I-(Δx/3)*(A+xs[i+1]*B)
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
    elseif step=="irk23"   #Two-stage third-order
        for i=1:xl-1
            rhs1=c0
            lhs1=I-Δx*(1/2+sqrt(3)/6)*(A+(xs[i]+Δx*(1/2+sqrt(3)/6))*B)
            temp1=lhs1\rhs1
            rhs2=c0-(sqrt(3)/3)*Δx*(A+(xs[i]+Δx*(1/2+sqrt(3)/6))*B)*temp1
            lhs2=I-Δx*(1/2+sqrt(3)/6)*(A+(xs[i]+Δx*(1/2-sqrt(3)/6))*B)
            temp2=lhs2\rhs2
            c0=c0+(Δx/2)*((A+(xs[i]+Δx*(1/2+sqrt(3)/6))*B)*temp1+(A+(xs[i]+Δx*(1/2-sqrt(3)/6))*B)*temp2)
            c0_pdf=(A+xs[i+1]*B)*c0
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1] = c0[end]
                TW_pdf[i+1] = c0_pdf[end]
            end
        end
    elseif step=="irk33"    #Three-stage third-order
        for i = 1:xl - 1
            rhs1 = c0 
            lhs1 = I-x33*Δx*(A+(xs[i]+x33*Δx)*B)
            temp1 = lhs1\rhs1
            rhs2 = c0+Δx*((1-x33)/2)*(A+(xs[i]+x33*Δx)*B)*temp1
            lhs2 = I-x33*Δx*(A+(xs[i]+((1+x33)/2)*Δx)*B)
            temp2 = lhs2\rhs2
            rhs3 = c0+Δx*(-3*x33^2/2+4*x33-1/4)*(A+(xs[i]+x33*Δx)*B)*temp1+Δx*(3*x33^2/2-5*x33+5/4)*(A+(xs[i]+((1+x33)/2)*Δx)*B)*temp2
            lhs3 = I-x33*Δx*(A+xs[i+1]*B)
            c0 = lhs3\rhs3
            c0_pdf=(A+xs[i+1]*B)*c0
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1] = c0[end]
                TW_pdf[i+1] = c0_pdf[end]
            end
        end
    elseif step=="irk43"     #Four-stage third-order
        for i=1:xl-1
            rhs1=c0
            lhs1=I-(Δx/2)*(A+(xs[i]+Δx/2)*B)
            temp1=lhs1\rhs1
            rhs2=c0+(Δx/6)*(A+(xs[i]+Δx/2)*B)*temp1
            lhs2=I-(Δx/2)*(A+(xs[i]+2*Δx/3)*B)
            temp2=lhs2\rhs2
            rhs3=c0-(Δx/2)*(A+(xs[i]+Δx/2)*B)*temp1+(Δx/2)*(A+(xs[i]+2*Δx/3)*B)*temp2
            lhs3=I-(Δx/2)*(A+(xs[i]+Δx/2)*B)
            temp3=lhs3\rhs3
            rhs4=c0+(3*Δx/2)*(A+(xs[i]+Δx/2)*B)*temp1-(3*Δx/2)*(A+(xs[i]+2*Δx/3)*B)*temp2+(Δx/2)*(A+(xs[i]+Δx/2)*B)*temp3
            lhs4=I-(Δx/2)*(A+xs[i+1]*B)
            c0=lhs4\rhs4
            c0_pdf=(A+xs[i+1]*B)*c0
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1] = c0[end]
                TW_pdf[i+1] = c0_pdf[end]
            end
        end
    elseif step=="irk34"       #Three-stage fourth-order
        for i=1:xl-1
            rhs1 = c0
            lhs1 = I-Δx*((1+α)/2)*(A+(xs[i]+((1+α)/2)*Δx)*B)
            temp1 = lhs1\rhs1
            rhs2 = c0-(α/2)*Δx*(A+(xs[i]+Δx*((1+α)/2))*B)*temp1
            lhs2 = I-Δx*((1+α)/2)*(A+(xs[i]+Δx/2)*B)
            temp2 = lhs2\rhs2
            rhs3 = c0+Δx*(1+α)*(A+(xs[i]+Δx*((1+α)/2))*B)*temp1-Δx*(1+2*α)*(A+(xs[i]+Δx/2)*B)*temp2
            lhs3 = I-Δx*((1+α)/2)*(A+(xs[i]+((1-α)/2)*Δx)*B)
            temp3 = lhs3\rhs3
            c0 = c0+(Δx/(6*α^2))*(A+(xs[i]+((1+α)/2)*Δx)*B)*temp1+Δx*(1-(1/(3*α^2)))*(A+(xs[i]+Δx/2)*B)*temp2+(Δx/(6*α^2))*(A+(xs[i]+((1-α)/2)*Δx)*B)*temp3
            c0_pdf=(A+xs[i+1]*B)*c0    
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1] = c0[end]
                TW_pdf[i+1] = c0_pdf[end]
            end
        end
    elseif step=="irkn34"       #Three-stage fourth-order 
        for i=1:xl-1 
            rhs1 = c0
            lhs1 = I-Δx*x34*(A+(xs[i]+x34*Δx)*B)
            temp1 = lhs1\rhs1
            rhs2 = c0+(1/2-x34)*Δx*(A+(xs[i]+Δx*x34)*B)*temp1
            lhs2 = I-Δx*x34*(A+(xs[i]+Δx/2)*B)
            temp2 = lhs2\rhs2
            rhs3 = c0+Δx*2*x34*(A+(xs[i]+Δx*x34)*B)*temp1+Δx*(1-4*x34)*(A+(xs[i]+Δx/2)*B)*temp2
            lhs3 = I-Δx*x34*(A+(xs[i]+(1-x34)*Δx)*B)
            temp3 = lhs3\rhs3
            c0 = c0+(Δx/(6*(1-2*x34)^2))*(A+(xs[i]+x34*Δx)*B)*temp1+Δx*(((3*(1-2*x34)^2)-1)/(3*(1-2*x34)^2))*(A+(xs[i]+Δx/2)*B)*temp2+(Δx/(6*(1-2*x34)^2))*(A+(xs[i]+(1-x34)*Δx)*B)*temp3
            c0_pdf = (A+xs[i+1]*B)*c0
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
    elseif step=="bdf4"
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
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1]=c0[end]
                TW_pdf[i+1]=c0_pdf[end]
            end
        end
    elseif step=="bdf5"
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
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1]=c0[end]
                TW_pdf[i+1]=c0_pdf[end]
            end
        end
    elseif step=="bdf6"
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
            if method=="spectral"
                TW_cdf[i+1] = sum(c0.*integ)|> real
                TW_pdf[i+1] = sum(c0_pdf.*integ)|> real
            elseif method=="finite"
                TW_cdf[i+1]=c0[end]
                TW_pdf[i+1]=c0_pdf[end]
            end
        end
    end
    return TW_cdf,TW_pdf
end
