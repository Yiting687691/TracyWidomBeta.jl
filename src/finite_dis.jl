function finite_dis(β,cheb_para;x0=13.0,xN=-10.0,Δx=-0.01,M=10^3)
    
    # Set up the domain and initial condition
    x=x0:Δx:xN;xl=length(x);
    θM=1;h=θM/(M-1);θ=0:h:θM;θ=θ[2:end];θ=θ*pi;h=h*pi;θM=θ[end];l=M-1;
    initial=ones(l,1);
    ind=convert(Int64,ceil((pi/2)/h));
    for i=1:ind-1
        initial[i]=cdf(Normal(),(x0-(cot(i*h))^2)/(sqrt((4/β)*cot(i*h))))[1]
    end
    final=initial;
    final_interest_cdf=zeros(xl,1);
    final_interest_pdf=zeros(xl,1);
    final_interest_cdf[1]=1;
    final_interest_pdf[1]=0;
    
    # Set up the discretization matrices
    p1=-2*ones(Int64,l-1,1);p1=vcat(p1,9/4);p2=ones(Int64,l-2,1);p2=vcat(p2,-6);p3=zeros(Int64,l-3,1);p3=vcat(p3,11/2);
    p4=zeros(Int64,l-4,1);p4=vcat(p4,-2);p5=zeros(Int64,l-5,1);p5=vcat(p5,1/4);
    T=spdiagm(0=>vec(p1),1=>fill(1.0,l-1),-1=>vec(p2),-2=>vec(p3),-3=>vec(p4),-4=>vec(p5));
    tt=(-2*(sin.(θ)).^4)/(β*(h^2));
    TT=spdiagm(0=>vec(tt))*T;
    e1=ones(Int64,l-2,1);e1=vcat(e1,4);e2=zeros(Int64,l-1,1);e2=vcat(e2,-3);e3=zeros(Int64,l-3,1);e3=vcat(e3,-1);
    U=spdiagm(0=>vec(e2),1=>fill(-1.0,l-1),-1=>vec(e1),-2=>vec(e3));
    
    # Step forward with trapezoidal rule
    for j=1:xl-1
        u1=(1/(2*h))*((x[j].+(2*sin.(2*θ))/β).*(sin.(θ)).^2-(cos.(θ)).^2);
        U1=spdiagm(0=>vec(u1))*U;
        u2=(1/(2*h))*((x[j+1].+(2*sin.(2*θ))/β).*(sin.(θ)).^2-(cos.(θ)).^2);
        U2=spdiagm(0=>vec(u2))*U;
        rhs=(I+(Δx/2)*TT+(Δx/2)*U1)*final;
        lhs=I-(Δx/2)*TT-(Δx/2)*U2;
        final=lhs\rhs;
        final_pdf=TT*final+U2*final;
        final_interest_pdf[j+1]=final_pdf[end];
        final_interest_cdf[j+1]=final[end];
    end
    
    # Interpolation
    ϕ = y -> (erf.(y) .+ 1.0)/2;
    j=PeriodicSegment(xN,x0);
    S = Laurent(j);
    final_int_pdf=final_interest_pdf[2:end] |> reverse;
    final_int_cdf=final_interest_cdf[2:end] |> reverse;
    xx=x[2:end] |> reverse;
    f_pdf = Fun(S, ApproxFun.transform(S,final_int_pdf))
    f_cdf = Fun(S, ApproxFun.transform(S,final_int_cdf - ϕ(xx)))
    S2 = x[end]..x[1] |> Chebyshev
    pdf_cheb = Fun(f_pdf,S2,cheb_para) |> real
    cdf_cheb = Fun(f_cdf,S2,cheb_para) + Fun(ϕ,S2) |> real
    return cdf_cheb,pdf_cheb
end
