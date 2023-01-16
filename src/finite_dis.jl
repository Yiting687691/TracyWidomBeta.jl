function fd(β;cheb,interp,pdf,x0,xN,Δx_f,M_f)

    # Set up the initial condition
    Φ=x -> erf.(x/sqrt(2))/2 + 0.5
    g=(x0,β,t) -> Φ((x0 - cot(t).^2)./sqrt.(4/β * cot.(t)))
    h0=(x0,β,t) -> t < pi/2 ? g(x0,β,t) : 1.0
    mgrid=(n,L) -> L*(1:n)/n
    θ=mgrid(M_f,pi)
    h=(1/M_f)*pi
    c0=map(t -> h0(x0,β,t),θ)
    
    # Set up the time domain
    xs = x0:Δx_f:xN |> Array
    if abs(xs[end])<abs(xN)
        xs=x0:Δx_f:(xN+Δx_f) |> Array
    end
    xl=length(xs);
    
    # Set up the discretization matrices
    T=spdiagm(0=>fill(-2.0,M_f),1=>fill(1.0,M_f-1),-1=>fill(1.0,M_f-1));
    tt=(-2*(sin.(θ)).^4)/(β*(h^2));
    TT=spdiagm(0=>vec(tt))*T;
    um1=ones(Int64,M_f-2,1);um1=vcat(um1,4);ud=zeros(Int64,M_f-1,1);ud=vcat(ud,-3);um2=zeros(Int64,M_f-3,1);um2=vcat(um2,-1);
    U=spdiagm(0=>vec(ud),1=>fill(-1.0,M_f-1),-1=>vec(um1),-2=>vec(um2));
    
    # Step forward with trapezoidal rule
    final=c0;
    final_interest_cdf=zeros(xl,1);
    final_interest_pdf=zeros(xl,1);
    final_interest_cdf[1]=1;
    final_interest_pdf[1]=0;
    for j=1:xl-1
        u1=(1/(2*h))*((xs[j].+(2*sin.(2*θ))/β).*(sin.(θ)).^2-(cos.(θ)).^2);
        U1=spdiagm(0=>vec(u1))*U;
        u2=(1/(2*h))*((xs[j+1].+(2*sin.(2*θ))/β).*(sin.(θ)).^2-(cos.(θ)).^2);
        U2=spdiagm(0=>vec(u2))*U;
        rhs=(I+(Δx_f/2)*TT+(Δx_f/2)*U1)*final;
        lhs=I-(Δx_f/2)*TT-(Δx_f/2)*U2;
        final=lhs\rhs;
        final_pdf=TT*final+U2*final;
        final_interest_pdf[j+1]=final_pdf[end];
        final_interest_cdf[j+1]=final[end];
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
