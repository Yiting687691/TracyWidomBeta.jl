function TW(β;cheb=10^3,method="finite",interp=true,pdf=false,x0=13.0,xN=-10.0,Δx_f=-0.01,Δx_s=-0.002,M_f=Int(floor(-10/(Δx_f))),M_s=Int(floor(-16/(Δx_s))),l=20)

    # Set up the time domain
    if method=="finite"
        xs = x0:Δx_f:xN |> Array
        if abs(xs[end])<abs(xN)
            xs=x0:Δx_f:(xN+Δx_f) |> Array
        end
    elseif method=="bdf3" || method=="bdf4" || method=="bdf5" || method=="bdf6"
        xs = x0:Δx_s:xN |> Array
        if abs(xs[end])<abs(xN)
            xs=x0:Δx_s:(xN+Δx_s) |> Array
        end
    end
    xl=length(xs)
    
    # Set up the initial condition
    if method=="finite"
        c0=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[1]
        h=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[2]
        θ=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[3]
    elseif method=="bdf3"
        c0=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[1]
        c1=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[2]
        c2=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[3]
        integ=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[4]
    elseif method=="bdf4"
        c0=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[1]
        c1=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[2]
        c2=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[3]
        c3=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[4]
        integ=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[5]
    elseif method=="bdf5"
        c0=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[1]
        c1=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[2]
        c2=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[3]
        c3=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[4]
        c4=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[5]
        integ=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[6]
    elseif method=="bdf6"
        c0=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[1]
        c1=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[2]
        c2=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[3]
        c3=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[4]
        c4=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[5]
        c5=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[6]
        integ=initial_gen(β;method,x0,Δx_s,M_f,M_s,l)[7]
    end

    # Set up the matrices
    A=matrix_gen(β;method,M_f,M_s,l)[1]
    B=matrix_gen(β;method,M_f,M_s,l)[2]

    # Step forward in time
    if method=="finite"
        final=c0
        final_interest_cdf=zeros(xl,1)
        final_interest_pdf=zeros(xl,1)
        final_interest_cdf[1]=1
        final_interest_pdf[1]=0
        tt=(-2*(sin.(θ)).^4)/(β*(h^2))
        TT=spdiagm(0=>vec(tt))*A
        for j=1:xl-1
            u1=(1/(2*h))*((xs[j].+(2*sin.(2*θ))/β).*(sin.(θ)).^2-(cos.(θ)).^2)
            U1=spdiagm(0=>vec(u1))*B
            u2=(1/(2*h))*((xs[j+1].+(2*sin.(2*θ))/β).*(sin.(θ)).^2-(cos.(θ)).^2)
            U2=spdiagm(0=>vec(u2))*B
            rhs=(I+(Δx_f/2)*TT+(Δx_f/2)*U1)*final
            lhs=I-(Δx_f/2)*TT-(Δx_f/2)*U2
            final=lhs\rhs
            final_pdf=TT*final+U2*final
            final_interest_pdf[j+1]=final_pdf[end]
            final_interest_cdf[j+1]=final[end]
        end
    elseif method=="bdf3"
        RR = (A+(x0*B))*(6*Δx_s)
        IRR = 11I - RR
        ΔRR = (6Δx_s^2)*B
        final_interest_cdf=zeros(xl,1)
        final_interest_pdf=zeros(xl,1)
        final_interest_pdf[1] = 0
        final_interest_cdf[1] = 1
        for i = 2:length(xs)
            temp = 18c0 - 9c1 + 2c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            final_interest_cdf[i] = real(sum(c0.*integ))
            final_interest_pdf[i] = real(sum((A+xs[i]*B)*c0.*integ))
        end
    elseif method=="bdf4"
        RR = (A+(x0*B))*(12Δx_s)
        IRR = 25I - RR
        ΔRR = (12Δx_s^2)*B
        final_interest_cdf=zeros(xl,1)
        final_interest_pdf=zeros(xl,1)
        final_interest_pdf[1] = 0
        final_interest_cdf[1] = 1
        for i = 2:length(xs)
            temp = 48c0 - 36c1 + 16c2 - 3c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            final_interest_cdf[i] = real(sum(c0.*integ))
            final_interest_pdf[i] = real(sum((A+xs[i]*B)*c0.*integ))
        end
    elseif method=="bdf5"
        RR = (A+(x0*B))*(60Δx_s)
        IRR = 137I - RR
        ΔRR = (60Δx_s^2)*B
        final_interest_cdf=zeros(xl,1)
        final_interest_pdf=zeros(xl,1)
        final_interest_pdf[1] = 0
        final_interest_cdf[1] = 1
        for i = 2:length(xs)
            temp = 300c0 - 300c1 + 200c2 - 75c3 + 12c4
            c4 = c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            final_interest_cdf[i] = real(sum(c0.*integ))
            final_interest_pdf[i] = real(sum((A+xs[i]*B)*c0.*integ))
        end
    elseif method=="bdf6"
        RR = (A+(x0*B))*(60Δx_s)
        IRR = 147I - RR
        ΔRR = (60Δx_s^2)*B
        final_interest_cdf=zeros(xl,1)
        final_interest_pdf=zeros(xl,1)
        final_interest_pdf[1] = 0
        final_interest_cdf[1] = 1
        for i = 2:length(xs)
            temp = 360c0 - 450c1 + 400c2 - 225c3 + 72c4 - 10c5
            c5 = c4
            c4 = c3
            c3 = c2
            c2 = c1
            c1 = c0
            IRR -= ΔRR
            c0 = IRR\temp
            final_interest_cdf[i] = real(sum(c0.*integ))
            final_interest_pdf[i] = real(sum((A+xs[i]*B)*c0.*integ))
        end
    end

    # Interpolation
    ϕ = y -> (erf.(y) .+ 1.0)/2;
    j=PeriodicSegment(xs[end],xs[1]);
    S = Laurent(j);
    final_int_pdf=final_interest_pdf[2:end] |> reverse;
    final_int_cdf=final_interest_cdf[2:end] |> reverse;
    xx=xs[2:end] |> reverse;
    f_pdf = Fun(S, ApproxFun.transform(S,final_int_pdf))
    f_cdf = Fun(S, ApproxFun.transform(S,final_int_cdf - ϕ(xx)))
    S2 = xs[end]..xs[1] |> Chebyshev
    pdf_cheb = Fun(f_pdf,S2,cheb) |> real
    cdf_cheb = Fun(f_cdf,S2,cheb) + Fun(ϕ,S2) |> real
    
    if interp && pdf==false
        return cdf_cheb
    elseif interp && pdf
        return pdf_cheb
    elseif interp==false && pdf==false
        return xs,final_interest_cdf
    elseif interp==false && pdf
        return xs,final_interest_pdf
    else
        return "Input valid arguments"
    end
end
