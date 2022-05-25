beta=2;
S=sol(beta);
delta_theta = 0.001*pi
final_theta=3*pi-delta_theta; domain=0:delta_theta:final_theta; domain=domain[2:end-1];
x=domain[1:end-998]
initial_time = 13.0
final_time = -10.0
delta_x = -0.01
tt=initial_time:delta_x:final_time;
δ= 0.6
m = 9
slope = 8
dx = delta_theta
dt = -1
sc = 5.0
up = 3
C = 1

X1 = x[1]:dx:x[end]
ts = 4:dt:-4
X2 = ts.-2.7
X2 = X2/X2[end]*(δ*m)
X3 = X2/slope

d = zeros(length(ts),length(x))

for i=1:length(ts)
    posi=findall(x->x==ts[i],tt)[1]
    r=S[posi,:]
    r=r[1:end-998]
    d[i,:]=r
end


x_range = [x[1],x[end]+ δ/slope*m]

p3 = plot()

for j = m:-1:1
    plot!(x .+ δ*(j-1)/slope, sc*d[j,:] .+ δ*(j-1) .+ C, grid = false, legend = false, lw = 2, color = :black, fillrange = δ*(j-1), fillcolor = :white, seriesalpha = 2, alpha = 0.9)
end
#xlabel!("s")
#ylabel!("β")


## draws the x axis
plot!(X1, 0*X1, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]) )

## draws the t axis
plot!(X3 .+ x[end].+0.2, X2.+2, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]) )

## labelling for the t axis
X2 = ts.-2.7
X2 = X2/X2[end]*(δ*m)
X3 = X2/slope
for j = 1:length(X3)
    annotate!(X3[j].+ x[end].+0.2, X2[j]+2.1, text("-"))
    annotate!(X3[j].+ x[end].+0.25, X2[j]+2, text(ts[j], :left,6))
end

## labelling for the x axis
X1=x[1]:1:x[end]
for j = 1:length(X1)
    annotate!(X1[j], 0, text("|",7))
    annotate!(X1[j], -.4, text(round(X1[j]),8,:below))
end

display(p3)
