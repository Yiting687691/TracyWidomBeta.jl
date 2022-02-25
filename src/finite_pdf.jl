function finite_pdf(beta,s;initial_time = 13.0, final_time = -10.0, delta_x = -0.01, delta_theta = 0.001*pi)
    time=initial_time:delta_x:final_time;
    final_theta=3*pi-delta_theta; domain=0:delta_theta:final_theta; domain=domain[2:end-1];
    initial=ones(length(domain),1);
    ind=convert(Int64,ceil((pi/2)/delta_theta));
    for j=1:ind-1
        initial[j]=cdf(Normal(),(initial_time-(cot(j*delta_theta))^2)/(sqrt((4/beta)*cot(j*delta_theta))))
    end
    final=initial;
    g1=zeros(length(initial),1); g2=g1;
    final_in=zeros(length(time),1);
    final_interest=zeros(length(time),1); final_interest[1]=0;
    e=ones(Int64,length(domain)-1,1);
    p=ones(Int64,length(domain),1);
    T=sparse(Array(1:length(domain)-1),Array(2:length(domain)),vec(e),length(domain),length(domain))+
        sparse(Array(2:length(domain)),Array(1:length(domain)-1),vec(e),length(domain),length(domain))+
        sparse(Array(1:length(domain)),Array(1:length(domain)),-2*vec(p),length(domain),length(domain));
    kk=(-2*(sin.(domain)).^4)/(beta*(delta_theta^2));
    TT=sparse(Array(1:length(domain)),Array(1:length(domain)),kk,length(domain),length(domain))*T;
    U=sparse(Array(1:length(domain)-1),Array(2:length(domain)),-1*vec(e),length(domain),length(domain))+
        sparse(Array(2:length(domain)),Array(1:length(domain)-1),vec(e),length(domain),length(domain));
    for k=1:length(time)-1
        u1=(1/(2*delta_theta))*((time[k].+(2*sin.(2*domain))/beta).*(sin.(domain)).^2-(cos.(domain)).^2);
        U1=sparse(Array(1:length(domain)),Array(1:length(domain)),u1,length(domain),length(domain))*U;
        u2=(1/(2*delta_theta))*((time[k+1].+(2*sin.(2*domain))/beta).*(sin.(domain)).^2-(cos.(domain)).^2);
        U2=sparse(Array(1:length(domain)),Array(1:length(domain)),u2,length(domain),length(domain))*U;
        g1[end]=((-2*(sin(final_theta))^4)/(beta*(delta_theta)^2))-(((time[k]+(2*sin(2*final_theta))/beta)*
            (sin(final_theta))^2-(cos(final_theta))^2)/(2*delta_theta));
        g2[end]=((-2*(sin(final_theta))^4)/(beta*(delta_theta)^2))-(((time[k+1]+(2*sin(2*final_theta))/beta)*
            (sin(final_theta))^2-(cos(final_theta))^2)/(2*delta_theta));
        rhs=(sparse(Array(1:length(domain)),Array(1:length(domain)),vec(p),length(domain),length(domain))+(delta_x/2)*
            TT+(delta_x/2)*U1)*final+(delta_x/2)*(g1+g2);
        lhs=sparse(Array(1:length(domain)),Array(1:length(domain)),vec(p),length(domain),length(domain))-(delta_x/2)*TT-(delta_x/2)*U2;
        final=lhs\rhs;
        final_in=TT*final+U2*final+g2;
        final_interest[k+1]=final_in[round(Int,pi/delta_theta)];
    end
    S = Chebyshev(final_time..initial_time);
    n=length(final_interest);
    while true
        if mod(n,10)==0
            break
        end
        n=n-1;
    end
    m=Int64(n/10);
    n=length(final_interest);
    V=zeros(n,m);
    for k = 1:m
           V[:,k] = Fun(S,[zeros(k-1);1]).(time)
    end
    f = Fun(S,V\vec(final_interest));
    return f(s)
end
