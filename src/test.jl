b=3;
initial_time = 13.0;
final_time = -10.0;
delta_x = -0.001;
delta_theta = 0.001*pi;
t=initial_time:delta_x:final_time;
final_theta=pi; dom=0:delta_theta:final_theta;dom=dom[2:end];
initial=ones(length(dom),1);
ind=convert(Int64,ceil((pi/2)/delta_theta));
for j=1:ind
    initial[j]=cdf(Normal(),(initial_time-(cot(j*delta_theta))^2)/(sqrt((4/b)*cot(j*delta_theta))))
end
final=initial;
final_interest=zeros(length(t),1); final_interest[1]=1;
e=ones(Int64,length(dom)-1,1);
eee=ones(Int64,length(dom)-1,1);eee[end]=4;
p=ones(Int64,length(dom),1);pp=p;pp[end]=0;ppp=zeros(Int64,length(dom),1);ppp[end]=-3;
q=zeros(Int64,length(dom)-2,1);q[end]=-1;
T=sparse(Array(1:length(dom)-1),Array(2:length(dom)),vec(e),length(dom),length(dom))+
    sparse(Array(2:length(dom)),Array(1:length(dom)-1),vec(e),length(dom),length(dom))+
    sparse(Array(1:length(dom)),Array(1:length(dom)),-2*vec(p),length(dom),length(dom));
kk=(-2*(sin.(dom)).^4)/(b*(delta_theta^2));
TT=sparse(Array(1:length(dom)),Array(1:length(dom)),kk,length(dom),length(dom))*T;
U=sparse(Array(1:length(dom)-1),Array(2:length(dom)),-1*vec(e),length(dom),length(dom))+
    sparse(Array(2:length(dom)),Array(1:length(dom)-1),vec(eee),length(dom),length(dom))+
    sparse(Array(1:length(dom)),Array(1:length(dom)),vec(ppp),length(dom),length(dom))+
    sparse(Array(3:length(dom)),Array(1:length(dom)-2),vec(q),length(dom),length(dom));
I=sparse(Array(1:length(dom)),Array(1:length(dom)),vec(p),length(dom),length(dom));
for k=1:length(t)-1
    u1=(1/(2*delta_theta))*((t[k].+(2*sin.(2*dom))/b).*(sin.(dom)).^2-(cos.(dom)).^2);
    U1=sparse(Array(1:length(dom)),Array(1:length(dom)),u1,length(dom),length(dom))*U;
    u2=(1/(2*delta_theta))*((t[k+1].+(2*sin.(2*dom))/b).*(sin.(dom)).^2-(cos.(dom)).^2);
    U2=sparse(Array(1:length(dom)),Array(1:length(dom)),u2,length(dom),length(dom))*U;
    rhs=(I+(delta_x/2)*TT+(delta_x/2)*U1)*final;
    lhs=I-(delta_x/2)*TT-(delta_x/2)*U2;
    final=lhs\rhs;
    final_interest[k+1]=final[end];
end
ind2=findall(x->x==5, t)[1]
ff=final_interest[ind2:end];
tt=t[ind2:end-1];
final_in=zeros(length(ff)-1,1)
final_in[1]=0
for i=1:length(ff)-2
    final_in[i+1]=(ff[i]-2*ff[i+1]+ff[i+2])/(2*delta_x)
end
