function finite_cdf(beta;initial_time = 13.0, final_time = -10.0, delta_x = -0.001,delta_theta = 0.001*pi)
    time=initial_time:delta_x:final_time;
    final_theta=2*pi; domain=0:delta_theta:final_theta;
    initial=ones(length(domain),1);
    ind=convert(Int64,ceil((pi/2)/delta_theta));
    for j=1:ind-1
        initial[j]=cdf(Normal(),(initial_time-(cot(j*delta_theta))^2)/(sqrt((4/beta)*cot(j*delta_theta))))
    end
    final=initial;
    final_interest=zeros(length(time),1); final_interest[1]=1;
    e=ones(Int64,length(domain)-1,1);ee=e;ee[end]=0;eee=e;eee[end]=4;
    p=ones(Int64,length(domain),1);pp=p;pp[end]=0;ppp=zeros(Int64,length(domain),1);ppp[end]=-3;
    q=zeros(Int64,length(domain)-2,1);q[end]=-1;
    T=sparse(Array(1:length(domain)-1),Array(2:length(domain)),vec(e),length(domain),length(domain))+
        sparse(Array(2:length(domain)),Array(1:length(domain)-1),vec(ee),length(domain),length(domain))+
        sparse(Array(1:length(domain)),Array(1:length(domain)),-2*vec(pp),length(domain),length(domain));
    kk=(-2*(sin.(domain)).^4)/(beta*(delta_theta^2));
    TT=sparse(Array(1:length(domain)),Array(1:length(domain)),kk,length(domain),length(domain))*T;
    U=sparse(Array(1:length(domain)-1),Array(2:length(domain)),-1*vec(e),length(domain),length(domain))+
        sparse(Array(2:length(domain)),Array(1:length(domain)-1),vec(eee),length(domain),length(domain))+
        sparse(Array(1:length(domain)),Array(1:length(domain)),vec(ppp),length(domain),length(domain))+
        sparse(Array(3:length(domain)),Array(1:length(domain)-2),vec(q),length(domain),length(domain));
    I=sparse(Array(1:length(domain)),Array(1:length(domain)),vec(p),length(domain),length(domain));
    for k=1:length(time)-1
        u1=(1/(2*delta_theta))*((time[k].+(2*sin.(2*domain))/beta).*(sin.(domain)).^2-(cos.(domain)).^2);
        U1=sparse(Array(1:length(domain)),Array(1:length(domain)),u1,length(domain),length(domain))*U;
        u2=(1/(2*delta_theta))*((time[k+1].+(2*sin.(2*domain))/beta).*(sin.(domain)).^2-(cos.(domain)).^2);
        U2=sparse(Array(1:length(domain)),Array(1:length(domain)),u2,length(domain),length(domain))*U;
        rhs=(I+(delta_x/2)*TT+(delta_x/2)*U1)*final;
        lhs=I-(delta_x/2)*TT-(delta_x/2)*U2;
        final=lhs\rhs;
        final_interest[k+1]=final[round(Int,pi/delta_theta)];
    end
    ϕ = x -> (erf.(x) .+ 1.0)/2;
    j=PeriodicSegment(final_time,initial_time);
    S = Laurent(j);
    final_int=final_interest[2:end] |> reverse;
    t=time[2:end] |> reverse;
    ff = Fun(S, ApproxFun.transform(S,final_int - ϕ(t)))
    S2 = -10..13 |> Chebyshev
    cdf_cheb = Fun(ff,S2,1000) + Fun(ϕ,S2) |> real
    return cdf_cheb,final_int,t
end
