beta=5
F=TW(beta)
f=F'
n=50
k=1000
ei=zeros(k,1)
for j=1:k
    N=rand(Normal(0,sqrt(2)), n)
    C1=vec(zeros(n-1,1))
    for i=1:n-1
        C1[i]=rand(Chi(beta*(n-i)), 1)[1]
    end
    H=(1/sqrt(beta))*SymTridiagonal(N,C1)
    Ht=n^(1/6)*(H-2*sqrt(n)*I);
    ei[j]=findmax(eigvals(Ht))[1]
end
s0=histogram(ei,normed=true,legend=false)
plot!(f,legend=false)
