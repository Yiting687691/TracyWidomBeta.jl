function model_test(beta,n)
    F=TW(beta)
    FF=TW(beta,method="BDF4")
    f=F'
    ff=FF'
    k=10^6
    ei=zeros(k,1)
    for j=1:k
        N=rand(Normal(0,sqrt(2)), n)
        C1=vec(zeros(n-1,1))
        for i=1:n-1
            C1[i]=rand(Chi(beta*(n-i)), 1)[1]
        end
        H=(1/sqrt(beta))*SymTridiagonal(N,C1)
        ei[j]=n^(1/6)*(eigmax(H)[1]-2*sqrt(n))
    end
    s=Plots.histogram(ei,normed=true,label="Normed Histogram",xticks=([-5;-5:0.5:1;1],[-5;-5:0.5:1;1]),dpi=500,size=(1200,800))
    Plots.plot!(f,xlims=[-5,1],label="Finite Difference",lw=3,dpi=500,size=(1200,800))
    Plots.plot!(ff,xlims=[-5,1],label="Spectral",lw=3,dpi=500,size=(1200,800))
    return s
end
