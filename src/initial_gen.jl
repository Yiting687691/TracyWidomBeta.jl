function initial_gen(β;method,x0,Δx_s,M_f,M_s,l)
    Φ = x -> erf.(x/sqrt(2))/2 + 0.5
    dϕ = x -> exp.(-x^2/2)*1/sqrt(2pi)
    g = (x0,β,t) -> Φ( (x0 - cot(t).^2)./sqrt.(4/β * cot.(t)) )
    h0 = (x0,β,t) -> t < pi/2 ? g(x0,β,t) : 1.0
    chain = (x0,β,t) -> (x0 + 3*cot.(t).^2).*csc.(t).*sec.(t)./sqrt.(16/β * cot.(t))
    dg = (x0,β,t) -> t ≈ 0.0 ? 0.0 : dϕ((x0 - cot(t).^2)./sqrt(4/β * cot.(t))).*chain(x0,β,t)
    dh0 = (x0,β,t) -> t < pi/2 ? dg(x0,β,t) : 0.0
    ft=(c,m,M,j)->(1/j)*trapz(m,vec(c.*exp.(-2*im*M*m/l)))
    idiffvec = m -> -floor(m/2):1:floor((m-1)/2)
    if method=="finite"
        mgrid=(n,L) -> L*(1:n)/n
        θ=mgrid(M_f,pi)
        h=(1/M_f)*pi
        c0=map(t -> h0(x0,β,t),θ)
        return c0,h,θ
    else
        mgrid = (n,L) -> L*(0:n)/n
        θ=mgrid(M_s,l*pi)
        N=idiffvec(M_s+1)
        integ=zeros(ComplexF64,length(N))
        for j=1:length(N)
            if N[j]==0
                integ[j]=pi
            else
                integ[j]=l/(2*im*N[j])*(exp(2*im*N[j]*pi/l)-1)
            end
        end
        c0=map(t -> dh0(x0,β,t),θ)
        c1=map(t -> dh0(x0-Δx_s,β,t),θ)
        c2=map(t -> dh0(x0-2*Δx_s,β,t),θ)
        c3=map(t -> dh0(x0-3*Δx_s,β,t),θ)
        c4=map(t -> dh0(x0-4*Δx_s,β,t),θ)
        c5=map(t -> dh0(x0-5*Δx_s,β,t),θ)
        c0t=zeros(ComplexF64,M_s+1)
        c1t=zeros(ComplexF64,M_s+1)
        c2t=zeros(ComplexF64,M_s+1)
        c3t=zeros(ComplexF64,M_s+1)
        c4t=zeros(ComplexF64,M_s+1)
        c5t=zeros(ComplexF64,M_s+1)
        for i=1:length(N)
            c0t[i]=ft(c0,θ,N[i],l*pi)
            c1t[i]=ft(c1,θ,N[i],l*pi)
            c2t[i]=ft(c2,θ,N[i],l*pi)
            c3t[i]=ft(c3,θ,N[i],l*pi)
            c4t[i]=ft(c4,θ,N[i],l*pi)
            c5t[i]=ft(c5,θ,N[i],l*pi)
        end
    end
    if method=="bdf3"
        return c0t,c1t,c2t,integ
    elseif method=="bdf4"
        return c0t,c1t,c2t,c3t,integ
    elseif method=="bdf5"
        return c0t,c1t,c2t,c3t,c4t,integ
    elseif method=="bdf6"
        return c0t,c1t,c2t,c3t,c4t,c5t,integ
    end
end
