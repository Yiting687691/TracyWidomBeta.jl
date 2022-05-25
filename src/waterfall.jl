Beta = 1:1:10
initial_time = 5
final_time = -6.0
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

δ= 0.5
m = 10
slope = 3
dx = 1.0
dt = 1
sc = 4.0
up = 4
C = 1

X1 = -6:dx:5
ts = Beta[1]:dt:Beta[end]
X2 = ts
X2 = X2/X2[end]*(δ*m)
X3 = X2/slope


x_range = [-6,5+ δ/slope*m]

p = plot()

for j = m:-1:1
    plot!(x .+ δ*(j-1)/slope, sc*d[j,:] .+ δ*(j-1) .+ C, grid = false, legend = false, lw = 2, color = :black, fillrange = δ*(j-1), fillcolor = :white, seriesalpha = 1, alpha = 0.9)
end

#xlabel!("s")
#ylabel!("β")


## draws the x axis
plot!(X1, 0*X1, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]) )

## draws the t axis
plot!(X3 .+ 5 , X2.+0.5, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]) )

## labelling for the t axis
for j = 1:length(X3)
    annotate!(X3[j].+ 5, X2[j]+.55, text("-"))
    annotate!(X3[j].+ 5 + .2, X2[j]+.5, text(ts[j], :left,8))
end

## labelling for the x axis
for j = 1:length(X1)
    annotate!(X1[j], 0, text("|",7))
    annotate!(X1[j], -.4, text(X1[j],8,:below))
end

display(p)
