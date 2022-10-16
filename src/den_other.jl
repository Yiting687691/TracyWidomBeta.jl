#This function gives the limiting density distribution of the largest n eigenvalues of the beta-Hermite ensemble.
function den_other(beta,n;initial_time = 13.0, final_time = -10.0, delta_x = -0.01,delta_theta = 0.001*pi)
    ss=Plots.plot()
    for i=1:n
        t=initial_time:delta_x:final_time;
        final_theta=i*pi;
        d=0:delta_theta:final_theta;
        d=d[2:end];
        initial=ones(length(d),1);
        ind=convert(Int64,ceil((pi/2)/delta_theta));
        for j=1:ind-1
            initial[j]=cdf(Normal(),[(initial_time-(cot(j*delta_theta))^2)/(sqrt((4/beta)*cot(j*delta_theta)))])[1]
        end
        final=initial;
        final_interest=zeros(length(t),1); final_interest[1]=0;
        e=ones(length(d)-1,1); e1=ones(length(d)-1,1); e1[end]=-6; e2=ones(length(d)-1,1); e2[end]=4;
        p=-2*ones(length(d),1); p[end]=9/4;
        p1=zeros(length(d),1); p1[end]=-3;
        p2=ones(length(d),1);
        q=zeros(length(d)-2,1); q[end]=11/2; q1=zeros(length(d)-2,1); q1[end]=-1;
        s=zeros(length(d)-3,1); s[end]=-2;
        r=zeros(length(d)-4,1); r[end]=1/4;
        T=sparse(Array(1:length(d)-1),Array(2:length(d)),vec(e),length(d),length(d))+
        sparse(Array(2:length(d)),Array(1:length(d)-1),vec(e1),length(d),length(d))+
        sparse(Array(1:length(d)),Array(1:length(d)),vec(p),length(d),length(d))+
        sparse(Array(3:length(d)),Array(1:length(d)-2),vec(q),length(d),length(d))+
        sparse(Array(4:length(d)),Array(1:length(d)-3),vec(s),length(d),length(d))+
        sparse(Array(5:length(d)),Array(1:length(d)-4),vec(r),length(d),length(d));
        kk=(-2*(sin.(d)).^4)/(beta*(delta_theta^2));
        TT=sparse(Array(1:length(d)),Array(1:length(d)),kk,length(d),length(d))*T;
        U=sparse(Array(1:length(d)-1),Array(2:length(d)),-1*vec(e),length(d),length(d))+
        sparse(Array(2:length(d)),Array(1:length(d)-1),vec(e2),length(d),length(d))+
        sparse(Array(1:length(d)),Array(1:length(d)),vec(p1),length(d),length(d))+
        sparse(Array(3:length(d)),Array(1:length(d)-2),vec(q1),length(d),length(d));
        I=sparse(Array(1:length(d)),Array(1:length(d)),vec(p2),length(d),length(d));
        for k=1:length(t)-1
            u1=(1/(2*delta_theta))*((t[k].+(2*sin.(2*d))/beta).*(sin.(d)).^2-(cos.(d)).^2);
            U1=sparse(Array(1:length(d)),Array(1:length(d)),u1,length(d),length(d))*U;
            u2=(1/(2*delta_theta))*((t[k+1].+(2*sin.(2*d))/beta).*(sin.(d)).^2-(cos.(d)).^2);
            U2=sparse(Array(1:length(d)),Array(1:length(d)),u2,length(d),length(d))*U;
            rhs=(I+(delta_x/2)*TT+(delta_x/2)*U1)*final;
            lhs=I-(delta_x/2)*TT-(delta_x/2)*U2;
            final=lhs\rhs;
            final_in=TT*final+U2*final;
            final_interest[k+1]=final_in[end];
        end
        ind2=findall(x->x==5, t)[1]
        ff=final_interest[ind2:end];
        tt=t[ind2:end];
        if i==1
            Plots.plot!(tt,ff,lw=3,dpi=500,size=(1200,800),label="Largest")
        elseif i==2
            Plots.plot!(tt,ff,lw=3,dpi=500,size=(1200,800),label="Second Largest")
        elseif i==3
            Plots.plot!(tt,ff,lw=3,dpi=500,size=(1200,800),label="Third Largest")
        end
    end
    return ss
end
