p4 = plot()
Beta = 1:0.5:4
for i = 1:length(Beta)
    beta = Beta[i]
    F=finite_cdf(beta)
    finite=F[1]
    ff=finite'
    if beta==1 || beta==2 || beta==4
        plot!(ff,lw = 3,legend=false,c=:red)
    else
        plot!(ff,lw = 0.5,legend=false,c=:black)
    end
end
xlabel!("s")
ylabel!("F(β,s)")
