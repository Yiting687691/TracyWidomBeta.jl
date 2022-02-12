function BDF4(beta,s;initial_time = 13.0, final_time = -10.0, delta_x = -0.0005, delta_theta = 0.0025*pi, l = 20)

#initial_time=13; final_time=-10; delta_x=-0.0005;
    time=initial_time:delta_x:final_time;
#delta_theta=0.0025*pi; l=20; f
    final_theta=l*pi; domain=0:delta_theta:final_theta+delta_theta;
    p=(-(length(domain)-1)/2):1:(length(domain)-1)/2;
    initial=zeros(length(domain),1);
    ind=convert(Int64,ceil((pi/2)/delta_theta))
    for j=1:ind-1
        initial[j]=((csc(j*delta_theta)*sec(j*delta_theta))*(3*(cot(j*delta_theta))^2+initial_time)/(4*sqrt(cot(j*delta_theta)/beta)))*
        pdf(Normal(),(initial_time-(cot(j*delta_theta))^2)/(sqrt((4/beta)*cot(j*delta_theta))))
    end
    initial_trans=zeros(ComplexF64,length(initial));
    for k=1:length(initial)
        initial_trans[k]=(1/(l*pi))*trapz(domain,vec(initial.*exp.(-2*im*p[k]*domain/l)));
    end
    final_interest=zeros(length(time),1); final_interest[1]=1;
    ehl=ones(Int64(length(domain)-(l/2)),1);
    el=ones(Int64(length(domain)-l),1);
    e3hl=ones(Int64(length(domain)-(3*l/2)),1);
    e2l=ones(Int64(length(domain)-2*l),1);
    e=ones(Int64(length(domain)),1);
    p=p*2*im/l;
    S_phl=sparse(Array(1:length(domain)-l/2),Array(((l/2)+1):length(domain)),vec(ehl),length(domain),length(domain));
    S_mhl=sparse(((l/2)+1):length(domain),1:length(domain)-l/2,vec(ehl),length(domain),length(domain));
    S_pl=sparse(1:length(domain)-l,(l+1):length(domain),vec(el),length(domain),length(domain));
    S_ml=sparse((l+1):length(domain),1:length(domain)-l,vec(el),length(domain),length(domain));
    S_p3hl=sparse(((3*l/2)+1):length(domain),1:length(domain)-3*l/2,vec(e3hl),length(domain),length(domain));
    S_m3hl=sparse(1:length(domain)-3*l/2,((3*l/2)+1):length(domain),vec(e3hl),length(domain),length(domain));
    S_p2l=sparse(1:length(domain)-2*l,(2*l+1):length(domain),vec(e2l),length(domain),length(domain));
    S_m2l=sparse((2*l+1):length(domain),1:length(domain)-2*l,vec(e2l),length(domain),length(domain));
    D1=sparse(1:length(domain),1:length(domain),vec(p),length(domain),length(domain));
    D2=D1^2;
    I=sparse(1:length(domain),1:length(domain),vec(e),length(domain),length(domain));
    A=(-2/beta)*((3/8)*I-(1/4)*S_ml-(1/4)*S_pl+(1/16)*S_p2l+(1/16)*S_m2l)*D2+
    ((1/2)*I+(1/4)*S_ml+(1/4)*S_pl-(2/beta)*((1/(2*im))*S_ml-(1/(2*im))*S_pl)*((1/2)*I-(1/4)*S_pl-(1/4)*S_ml)-
    (8/beta)*((-1/(8*im))*S_p3hl+(1/(8*im))*S_m3hl+(3/(8*im))*S_mhl-(3/(8*im))*S_phl)*((1/2)*S_mhl+(1/2)*S_phl))*D1-
    (4/beta)*((1/2)*I-(1/4)*S_pl-(1/4)*S_ml)*((1/2)*S_pl+(1/2)*S_ml)-
    (4/beta)*((1/(2*im))*S_mhl-(1/(2*im))*S_phl)*((1/2)*S_phl+(1/2)*S_mhl)*((1/(2*im))*S_ml-(1/(2*im))*S_pl)-
    2*((1/(2*im))*S_mhl-(1/(2*im))*S_phl)*((1/2)*S_phl+(1/2)*S_mhl);
    B=((-1/2)*I+(1/4)*S_ml+(1/4)*S_pl)*D1-2*((1/(2*im))*S_mhl-(1/(2*im))*S_phl)*((1/2)*S_phl+(1/2)*S_mhl);
    final=initial_trans

    N=(-(length(domain)-1)/2):1:(length(domain)-1)/2;
    integ=zeros(ComplexF64,length(N));
    for j=1:length(N)
        if N[j]==0
            integ[j]=pi;
        else
            integ[j]=l/(2*im*N[j])*(exp(2*im*N[j]*pi/l)-1);
        end
    end
    delta_x=delta_x/100;
    final1=zeros(ComplexF64,length(final));
    final2=zeros(ComplexF64,length(final));
    final3=zeros(ComplexF64,length(final));
    for i=1:300
        Y1=final;
        Y2=final+(delta_x/2)*(A+time[1]*B)*Y1;
        Y3=final+(delta_x/2)*(A+(time[1]+delta_x/2)*B)*Y2;
        Y4=final+delta_x*(A+(time[1]+delta_x/2)*B)*Y3;
        final=final+(delta_x/6)*((A+time[1]*B)*Y1+2*(A+(time[1]+delta_x/2)*B)*Y2+2*(A+(time[1]+delta_x/2)*B)*Y3+
        (A+(time[1]+delta_x)*B)*Y4);
        if i==100
            final1=final;
        elseif i==200
            final2=final;
        elseif i==300
            final3=final;
        end
    end
    final_interest[2]=real(sum(final1.*integ));
    final_interest[3]=real(sum(final2.*integ));
    final_interest[4]=real(sum(final3.*integ));
    delta_x=delta_x*100;
    finals = hcat(final,final1,final2,final3)
    # for k=5:length(time)
    #     #lhs=25*I-12*delta_x*(A+time[k]*B);
    #     lhs=25*I- (12*delta_x)*A - (12*delta_x*time[k])*B;
    #     rhs=48*final3-36*final2+16*final1-3*final;
    #     final4=lhs\rhs;
    #     final_interest[k]=real(sum(final4.*integ));
    #     final=final1;
    #     final1=final2;
    #     final2=final3;
    #     final3=final4;
    # end
    for k=5:length(time)
        final_interest[k]= one_step!(finals,delta_x,time[k],A,B,integ)
    end
    #print(final_interest[findall(x->x==s, time)][1])
    return final_interest, A, B
end

function one_step!(final::Array{ComplexF64},delta_x::Float64,timek::Float64,A::SparseMatrixCSC,B::SparseMatrixCSC,integ)
    lhs=25*I- (12*delta_x)*A - (12*delta_x*timek)*B;
    rhs=48*final[:,4]-36*final[:,3]+16*final[:,2]-3*final[:,1];
    final4=lhs\rhs;
    fi =real(sum(final4.*integ));
    final[:,1]=final[:,2];
    final[:,2]=final[:,3];
    final[:,3]=final[:,4];
    final[:,4]=final4;
    fi
end
