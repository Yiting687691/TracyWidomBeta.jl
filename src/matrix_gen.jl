function matrix_gen(β;method,M_f,M_s,l)
    diffvec = (L,m,j) -> ((-floor(m/2):1:floor((m-1)/2))*(2*1im*pi/L)).^j
    function 𝒟(L,m,j)
        spdiagm(diffvec(L,m,j))
    end
    if method=="finite"
        T=spdiagm(0=>fill(-2.0,M_f),1=>fill(1.0,M_f-1),-1=>fill(1.0,M_f-1))
        um1=ones(Int64,M_f-2,1);um1=vcat(um1,4);ud=zeros(Int64,M_f-1,1);ud=vcat(ud,-3);um2=zeros(Int64,M_f-3,1);um2=vcat(um2,-1)
        U=spdiagm(0=>vec(ud),1=>fill(-1.0,M_f-1),-1=>vec(um1),-2=>vec(um2))
        return T,U
    elseif method=="bdf3" || method=="bdf4" || method=="bdf5" || method=="bdf6"
        ll=convert(Int64,l/2)
        mme = spdiagm( ll => fill(1.0,M_s+1-ll), ll-M_s-1 => fill(1.0,ll))
        me = spdiagm( -ll => fill(1.0,M_s+1-ll), M_s+1-ll => fill(1.0,ll))
        ms = (me - mme)/2im;
        ms2 = (me^2 - mme^2)/2im;
        mc = (me + mme)/2;
        mc2 = (me^2 + mme^2)/2;
        DD = 𝒟(l*pi,M_s+1,1) |> sparse;
        A = (-2/β)*ms^4*DD^2 - (8/β)*ms^3*mc*DD - (2/β)*ms2*ms^2*DD + mc^2*DD - (4/β)*(ms*mc*ms2 + ms^2*mc2) - 2*mc*ms;
        B = -ms^2*DD - 2*ms*mc
        return A,B
    end
end
