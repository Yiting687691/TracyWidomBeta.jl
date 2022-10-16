#This function will return a matrix, where each row represents the density of H with respect to theta. The first row represents the initial density.
function density_theta(beta,final_theta;initial_time = 13.0, final_time = -10.0, delta_x = -0.001,delta_theta = 0.001*pi)
    t=initial_time:delta_x:final_time;
    d=0:delta_theta:final_theta;
    d=d[2:end];
    initial=zeros(length(d),1);
    ind=convert(Int64,ceil((pi/2)/delta_theta));
    for j=1:ind
        initial[j]=((csc(j*delta_theta)*sec(j*delta_theta))*(3*(cot(j*delta_theta))^2+initial_time)/(4*sqrt(cot(j*delta_theta)/beta)))*
        pdf(Normal(),(initial_time-(cot(j*delta_theta))^2)/(sqrt((4/beta)*cot(j*delta_theta))))
    end
    final=initial;
    final_interest=zeros(Float64,length(t),length(initial));
    final_interest[1,:]=initial';
    e=ones(Int64,length(d)-1,1);
    p=ones(Int64,length(d),1);
    eee=ones(Int64,length(d)-1,1);eee[end]=4;
    ppp=zeros(Int64,length(d),1);ppp[end]=-3;
    q=zeros(Int64,length(d)-2,1);q[end]=-1;
    T=sparse(Array(1:length(d)-1),Array(2:length(d)),vec(e),length(d),length(d))+
        sparse(Array(2:length(d)),Array(1:length(d)-1),vec(e),length(d),length(d))+
        sparse(Array(1:length(d)),Array(1:length(d)),-2*vec(p),length(d),length(d));
    kk=(-2*(sin.(d)).^4)/(beta*(delta_theta^2));
    TT=sparse(Array(1:length(d)),Array(1:length(d)),kk,length(d),length(d))*T;
    U=sparse(Array(1:length(d)-1),Array(2:length(d)),-1*vec(e),length(d),length(d))+
        sparse(Array(2:length(d)),Array(1:length(d)-1),vec(eee),length(d),length(d))+
        sparse(Array(1:length(d)),Array(1:length(d)),vec(ppp),length(d),length(d))+
        sparse(Array(3:length(d)),Array(1:length(d)-2),vec(q),length(d),length(d));
    I=sparse(Array(1:length(d)),Array(1:length(d)),vec(p),length(d),length(d));
    for k=1:length(t)-1
        u1=(1/(2*delta_theta))*(8*((sin.(d)).^3).*(cos.(d))/beta+(t[k].+(2*sin.(2*d))/beta).*(sin.(d)).^2-(cos.(d)).^2);
        U1=sparse(Array(1:length(d)),Array(1:length(d)),u1,length(d),length(d))*U;
        u2=(1/(2*delta_theta))*(8*((sin.(d)).^3).*(cos.(d))/beta+(t[k+1].+(2*sin.(2*d))/beta).*(sin.(d)).^2-(cos.(d)).^2);
        U2=sparse(Array(1:length(d)),Array(1:length(d)),u2,length(d),length(d))*U;
        i1=-(2*t[k]+1)*sin.(d).*cos.(d)-(2/beta)*((2*(sin.(d)).^2).*cos.(2*d)+2*sin.(d).*cos.(d).*sin.(2*d));
        I1=sparse(Array(1:length(d)),Array(1:length(d)),i1,length(d),length(d))*I;
        i2=-(2*t[k+1]+1)*sin.(d).*cos.(d)-(2/beta)*((2*(sin.(d)).^2).*cos.(2*d)+2*sin.(d).*cos.(d).*sin.(2*d));
        I2=sparse(Array(1:length(d)),Array(1:length(d)),i2,length(d),length(d))*I;
        rhs=(I+(delta_x/2)*TT+(delta_x/2)*U1+(delta_x/2)*I1)*final;
        lhs=I-(delta_x/2)*TT-(delta_x/2)*U2-(delta_x/2)*I2;
        final=lhs\rhs;
        final_interest[k+1,:]=final';
    end
    return final_interest,d,t
end
