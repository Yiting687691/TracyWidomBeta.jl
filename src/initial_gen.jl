function initial_gen(β;method,x0,Δx_f,Δx_s,M_f,M_s,l)
    Φ = x -> erf.(x/sqrt(2))/2 + 0.5
    dϕ = x -> exp.(-x^2/2)*1/sqrt(2pi)
    g = (x0,β,t) -> Φ( (x0 - cot(t).^2)./sqrt.(4/β * cot.(t)) )
    h0 = (x0,β,t) -> t < pi/2 ? g(x0,β,t) : 1.0
    chain = (x0,β,t) -> (x0 + 3*cot.(t).^2).*csc.(t).*sec.(t)./sqrt.(16/β * cot.(t))
    dg = (x0,β,t) -> t ≈ 0.0 ? 0.0 : dϕ((x0 - cot(t).^2)./sqrt(4/β * cot.(t))).*chain(x0,β,t)
    dh0 = (x0,β,t) -> t < pi/2 ? dg(x0,β,t) : 0.0
    mfft = x -> fftshift(fft(x))/length(x)

    if method=="finite"
        mgrid=(n,L) -> L*(1:n)/n
        θ=mgrid(M_f,pi)
        h=(1/M_f)*pi
        c0=map(t -> h0(x0,β,t),θ)
        c1=map(t -> h0(x0-Δx_f,β,t),θ)
        c2=map(t -> h0(x0-2*Δx_f,β,t),θ)
        c3=map(t -> h0(x0-3*Δx_f,β,t),θ)
        c4=map(t -> h0(x0-4*Δx_f,β,t),θ)
        c5=map(t -> h0(x0-5*Δx_f,β,t),θ)
    elseif method=="spectral"
        mgrid = (n,L) -> 2*L*(0:n-1)/n
        θ=mgrid(M_s,l*pi)
        h=(1/M_s)*2*l*pi
        c0=map(t -> dh0(x0,β,t),θ)|>mfft
        c1=map(t -> dh0(x0-Δx_s,β,t),θ)|>mfft
        c2=map(t -> dh0(x0-2*Δx_s,β,t),θ)|>mfft
        c3=map(t -> dh0(x0-3*Δx_s,β,t),θ)|>mfft
        c4=map(t -> dh0(x0-4*Δx_s,β,t),θ)|>mfft
        c5=map(t -> dh0(x0-5*Δx_s,β,t),θ)|>mfft
    end
    return c0,c1,c2,c3,c4,c5,h,θ
end
