p4 = plot()
Beta = 1:0.1:2
for i = 1:length(Beta)
    beta = Beta[i]
    F=finite_cdf(beta)
    finite=F[1]
    ff=finite'
    plot!(ff,lw = 1,legend=false)
end
xlabel!("s")
ylabel!("F(β,s)")
