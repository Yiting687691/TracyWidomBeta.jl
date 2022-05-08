Beta = 1:1:10
initial_time = 13.0
final_time = -10.0
delta_x = -0.01
x=initial_time:delta_x:final_time
d = zeros(length(Beta),length(x))
for i = 1:length(Beta)
    beta = Beta[i]
    F=finite_cdf(beta)
    finite=F[1]
    ff=finite'
    for j=1:length(x)
        t=x[j]
        d[i,j] = ff(t);
    end
end

δ= 1
m = 10
slope = 2.5
dx = 1.0
dt = 0.1
sc = 4.0
up = 4
C = 1

X1 = x[1]:dx:x[end]
ts = Beta[1]:dt:Beta[end]
X2 = ts .- Beta[1]
X2 = X2/X2[end]*(δ*m)
X3 = X2/slope


x_range = [x[1],x[end] + δ/slope*m]

p = plot()

for j = m:-1:1
    plot!(x .+ δ*(j-1)/slope, sc*d[j,:] .+ δ*(j-1) .+ C, grid = false, legend = false, lw = 2, color = :black, fillrange = δ*(j-1), fillcolor = :white, seriesalpha = 2, alpha = 1.9)
end
xlabel!("s")
ylabel!("β")

p

## draws the x axis
plot!(X1, 0*X1, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]) )

## draws the t axis
plot!(X3 .+ x[end] , X2, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]) )

## labelling for the t axis
for j = 1:length(X3)
    annotate!(X3[j].+ x[end], X2[j]+.1, text("-"))
    annotate!(X3[j].+ x[end] + .08, X2[j], text(ts[j], :left,8))
end

## labelling for the x axis
for j = 1:length(X1)
    annotate!(X1[j], 0, text("|",7))
    annotate!(X1[j], -.4, text(X1[j],8,:below))
end

display(p)
