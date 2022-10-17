function contour_H(beta;initial_time = 13.0, final_time = -10.0, delta_x = -0.001,delta_theta = 0.001*pi)
    time=initial_time:delta_x:final_time;
    final_theta=pi; domain=0:delta_theta:final_theta;domain=domain[2:end];
    initial=ones(length(domain),1);
    ind=convert(Int64,ceil((pi/2)/delta_theta));
    for j=1:ind
        initial[j]=cdf(Normal(),(initial_time-(cot(j*delta_theta))^2)/(sqrt((4/beta)*cot(j*delta_theta))))
    end
    final=initial;
    F=zeros(length(time),length(final));
    F[end,:]=final';
    final_interest=zeros(length(time),1); final_interest[1]=1;
    e=ones(Int64,length(domain)-1,1);
    f=ones(Int64,length(domain)-1,1);
    f[end]=4;
    p=ones(Int64,length(domain),1);
    ppp=zeros(Int64,length(domain),1);ppp[end]=-3;
    q=zeros(Int64,length(domain)-2,1);q[end]=-1;
    T=sparse(Array(1:length(domain)-1),Array(2:length(domain)),vec(e),length(domain),length(domain))+
        sparse(Array(2:length(domain)),Array(1:length(domain)-1),vec(e),length(domain),length(domain))+
        sparse(Array(1:length(domain)),Array(1:length(domain)),-2*vec(p),length(domain),length(domain));
    kk=(-2*(sin.(domain)).^4)/(beta*(delta_theta^2));
    TT=sparse(Array(1:length(domain)),Array(1:length(domain)),kk,length(domain),length(domain))*T;
    U=sparse(Array(1:length(domain)-1),Array(2:length(domain)),-1*vec(e),length(domain),length(domain))+
        sparse(Array(2:length(domain)),Array(1:length(domain)-1),vec(f),length(domain),length(domain))+
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
        F[end-k,:]=final';
        final_interest[k+1]=final[end];
    end
    t=reverse(time);
    p = Plots.contour(domain, t, F, fill=true,dpi=1000)
    return p
end
