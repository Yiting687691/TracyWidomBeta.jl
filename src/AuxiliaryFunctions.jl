function time_gen(x1,x2,Δx)
    xs=x1:Δx:x2 |> Array
    if abs(xs[end])<abs(x2)
        xs=x1:Δx:(x2+Δx) |> Array
    end
    return xs
end

function initial_gen(β;method,x0,Δx_f,Δx_s,M_f,M_s,l)
    Φ = x -> erf.(x/sqrt(2))/2 + 0.5
    dϕ = x -> exp.(-x^2/2)*1/sqrt(2pi)
    g = (x0,β,t) -> Φ( (x0 - cot(t).^2)./sqrt.(4/β * cot.(t)) )
    h0 = (x0,β,t) -> t < pi/2 ? g(x0,β,t) : 1.0
    chain = (x0,β,t) -> (x0 + 3*cot.(t).^2).*csc.(t).*sec.(t)./sqrt.(16/β * cot.(t))
    dg = (x0,β,t) -> t ≈ 0.0 ? 0.0 : dϕ((x0 - cot(t).^2)./sqrt(4/β * cot.(t))).*chain(x0,β,t)
    dh0 = (x0,β,t) -> t < pi/2 ? dg(x0,β,t) : 0.0
    mfft = x -> fftshift(fft(x))/length(x)

    if method=="finite"
        mgrid=(n,L) -> L*(1:n)/n
        θ=mgrid(M_f,pi)
        h=(1/M_f)*pi
        c0=map(t -> h0(x0,β,t),θ)
        c1=map(t -> h0(x0-Δx_f,β,t),θ)
        c2=map(t -> h0(x0-2*Δx_f,β,t),θ)
        c3=map(t -> h0(x0-3*Δx_f,β,t),θ)
        c4=map(t -> h0(x0-4*Δx_f,β,t),θ)
        c5=map(t -> h0(x0-5*Δx_f,β,t),θ)
    elseif method=="spectral"
        mgrid = (n,L) -> 2*L*(0:n-1)/n
        θ=mgrid(M_s,l*pi)
        h=(1/M_s)*2*l*pi
        c0=map(t -> dh0(x0,β,t),θ)|>mfft
        c1=map(t -> dh0(x0-Δx_s,β,t),θ)|>mfft
        c2=map(t -> dh0(x0-2*Δx_s,β,t),θ)|>mfft
        c3=map(t -> dh0(x0-3*Δx_s,β,t),θ)|>mfft
        c4=map(t -> dh0(x0-4*Δx_s,β,t),θ)|>mfft
        c5=map(t -> dh0(x0-5*Δx_s,β,t),θ)|>mfft
    end
    return c0,c1,c2,c3,c4,c5,h,θ
end

function matrix_gen(β;method,M_f,M_s,h,θ,l)
    diffvec = (L,m,j) -> ((-floor(m/2):1:floor((m-1)/2))*(1im*pi/L)).^j
    function 𝒟(L,m,j)
        spdiagm(diffvec(L,m,j))
    end
    if method=="finite"
        T=spdiagm(0=>fill(-2.0,M_f),1=>fill(1.0,M_f-1),-1=>fill(1.0,M_f-1))
        tt=(-2*(sin.(θ)).^4)/(β*(h^2))
        T=spdiagm(0=>vec(tt))*T
        um1=ones(Int64,M_f-2,1);um1=vcat(um1,4);ud=zeros(Int64,M_f-1,1);ud=vcat(ud,-3);um2=zeros(Int64,M_f-3,1);um2=vcat(um2,-1)
        U=spdiagm(0=>vec(ud),1=>fill(-1.0,M_f-1),-1=>vec(um1),-2=>vec(um2))
        tt2=(1/(2*h))*((2*sin.(2*θ)/β).*(sin.(θ)).^2-(cos.(θ)).^2)
        A=T+spdiagm(0=>vec(tt2))*U
        uu=(1/(2*h))*(sin.(θ)).^2
        B=spdiagm(0=>vec(uu))*U
    elseif method=="spectral"
        mme = spdiagm( l => fill(1.0,M_s-l), l-M_s => fill(1.0,l))
        me = spdiagm( -l => fill(1.0,M_s-l), M_s-l => fill(1.0,l))
        ms = (me - mme)/2im;
        ms2 = (me^2 - mme^2)/2im;
        mc = (me + mme)/2;
        mc2 = (me^2 + mme^2)/2;
        DD = 𝒟(l*pi,M_s,1) |> sparse;
        A = (-2/β)*ms^4*DD^2 - (8/β)*ms^3*mc*DD - (2/β)*ms2*ms^2*DD + mc^2*DD - (4/β)*(ms*mc*ms2 + ms^2*mc2) - 2*mc*ms;
        B = -ms^2*DD - 2*ms*mc
    end
    return A,B
end

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

function Fourier_interp(xs,cheb,TW_cdf,TW_pdf)  
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
