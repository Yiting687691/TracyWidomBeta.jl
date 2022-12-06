function spectral_dis(β,cheb_para;x0=13.0,xN=-10.0,Δx=-0.002,M=8000,l=20,period_flag=true)
    x=x0:Δx:xN;xl=length(x)
    θM=l;h=θM/(M-1);
    θ=0:h:θM;
    θ=θ*pi;
    MM=-floor(M/2):1:floor((M-1)/2)
    initial=zeros(M,1)
    initial1=zeros(M,1)
    initial2=zeros(M,1)
    initial3=zeros(M,1)
    initial4=zeros(M,1)
    ind=convert(Int64,ceil((pi/2)/h))
    for j=1:ind-1
        initial[j]=((csc(j*h)*sec(j*h))*(3*(cot(j*h))^2+x0)/(4*sqrt(cot(j*h)/β)))*
        pdf(Normal(),(x0-(cot(j*h))^2)/(sqrt((4/β)*cot(j*h))))
        initial1[j]=((csc(j*h)*sec(j*h))*(3*(cot(j*h))^2+x0-Δx)/(4*sqrt(cot(j*h)/β)))*
        pdf(Normal(),(x0-Δx-(cot(j*h))^2)/(sqrt((4/β)*cot(j*h))))
        initial2[j]=((csc(j*h)*sec(j*h))*(3*(cot(j*h))^2+x0-2*Δx)/(4*sqrt(cot(j*h)/β)))*
        pdf(Normal(),(x0-2*Δx-(cot(j*h))^2)/(sqrt((4/β)*cot(j*h))))
        initial3[j]=((csc(j*h)*sec(j*h))*(3*(cot(j*h))^2+x0-3*Δx)/(4*sqrt(cot(j*h)/β)))*
        pdf(Normal(),(x0-3*Δx-(cot(j*h))^2)/(sqrt((4/β)*cot(j*h))))
        initial4[j]=((csc(j*h)*sec(j*h))*(3*(cot(j*h))^2+x0-4*Δx)/(4*sqrt(cot(j*h)/β)))*
        pdf(Normal(),(x0-4*Δx-(cot(j*h))^2)/(sqrt((4/β)*cot(j*h))))
    end
    initial_trans=zeros(ComplexF64,M)
    initial_trans1=zeros(ComplexF64,M)
    initial_trans2=zeros(ComplexF64,M)
    initial_trans3=zeros(ComplexF64,M)
    initial_trans4=zeros(ComplexF64,M)
    for k=1:length(initial)
        initial_trans[k]=(1/θM)*trapz(θ,vec(initial.*exp.(-2*im*MM[k]*θ/l)))
        initial_trans1[k]=(1/θM)*trapz(θ,vec(initial1.*exp.(-2*im*MM[k]*θ/l)))
        initial_trans2[k]=(1/θM)*trapz(θ,vec(initial2.*exp.(-2*im*MM[k]*θ/l)))
        initial_trans3[k]=(1/θM)*trapz(θ,vec(initial3.*exp.(-2*im*MM[k]*θ/l)))
        initial_trans4[k]=(1/θM)*trapz(θ,vec(initial4.*exp.(-2*im*MM[k]*θ/l)))
    end
    final_interest_cdf=zeros(xl,1); final_interest_cdf[1]=1
    final_interest_pdf=zeros(xl,1); final_interest_pdf[1]=0
    ll=convert(Int64,l/2)
    if period_flag
        mme = spdiagm( ll => fill(1.0,M-ll), ll-M => fill(1.0,ll))
        me = spdiagm( -ll => fill(1.0,M-ll), M-ll => fill(1.0,ll))
    else
        mme = spdiagm( ll => fill(1.0,MM-ll))
        me = spdiagm( -ll => fill(1.0,MM-ll))
    end
    ms = (me - mme)/2im;
    ms2 = (me^2 - mme^2)/2im;
    mc = (me + mme)/2;
    mc2 = (me^2 + mme^2)/2;
    DD = 𝒟(θM,M,1) |> sparse;
    A = (-2/β)*ms^4*DD^2 - (8/β)*ms^3*mc*DD - (2/β)*ms2*ms^2*DD + mc^2*DD - (4/β)*(ms*mc*ms2 + ms^2*mc2) - 2*mc*ms;
    B = -ms^2*DD - 2*ms*mc;
    final=initial_trans
    final1=initial_trans1
    final2=initial_trans2
    final3=initial_trans3
    final4=initial_trans4
    integ=zeros(ComplexF64,length(MM))
    for j=1:length(MM)
        if MM[j]==0
            integ[j]=pi
        else
            integ[j]=l/(2*im*MM[j])*(exp(2*im*MM[j]*pi/l)-1)
        end
    end
    finals_cdf = hcat(final4,final3,final2,final1,final)
    finals_pdf = hcat(final4,final3,final2,final1,final)
    A1 = 137*I- (60*Δx)*A
    A2 = -(60*Δx)*B
    for k=2:xl
        final_interest_cdf[k]= one_step5!(finals_cdf,Δx,x[k],A1,A2,integ)
        final_interest_pdf[k]= one_step5_pdf!(finals_pdf,Δx,x[k],A1,A2,integ)
    end
    ϕ=y->(erf.(y).+1.0)/2
    j=PeriodicSegment(xN,x0)
    S = Laurent(j)
    final_int_cdf=final_interest_cdf[2:end] |> reverse
    final_int_pdf=final_interest_pdf[2:end] |> reverse
    x2=x[2:end] |> reverse
    f_cdf = Fun(S, ApproxFun.transform(S,final_int_cdf - ϕ(x2)))
    f_pdf = Fun(S, ApproxFun.transform(S,final_int_pdf))
    S2 = x[end]..x[1] |> Chebyshev
    cdf_cheb = Fun(f_cdf,S2,cheb_para) + Fun(ϕ,S2) |> real
    pdf_cheb = Fun(f_pdf,S2,cheb_para) |> real
    return cdf_cheb,pdf_cheb
end