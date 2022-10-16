#This function gives a waterfall plot of the solution with respect to theta from theta=0 to 10pi.
function water_den(beta,final_theta)
    result=density_theta(beta,final_theta)
    final_interest=result[1]
    d=result[2]
    t=result[3]
    x=d

    δ= 0.8
    m = 8
    slope = 2
    dt = -2
    sc = 200.0
    up = 10
    C = 1
    delt=(x[end]-x[1])/6
    X1=x[1]:delt:x[end]
    ts = 4:dt:-10
    X2 = ts.-4
    tick_shift = 0.01
    axis_shift = 0.97
    X2 = X2/X2[end]*(δ*(m-1)).+ tick_shift
    X3 = X2/slope
    dd = zeros(length(ts),length(x))
    for i=1:length(ts)
        posi=findall(x->x==ts[i],t)[1]
        r=final_interest[posi,:]
        dd[i,:]=r
    end
    x_range = [x[1],x[end]+ δ/slope*m]

    p=Plots.plot()
    for j = m:-1:1
        Plots.plot!(x .+ δ*(j-1)/slope, sc*dd[j,:] .+ δ*(j-1) .+ C, grid = false, legend = false, lw = 2, color = :black, fillrange = δ*(j-1), fillcolor = :white, seriesalpha = 2, alpha = 0.9,dpi=500,size=(1200,800))
    end

## draws the x axis
    Plots.plot!(X1, 0*X1, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]),dpi=500,size=(1200,800) )

## draws the t axis
    Plots.plot!(X3 .+ x[end].+0.1, X2.+axis_shift, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]),dpi=500,size=(1200,800) )

## labelling for the t axis
    for j = 1:length(X3)
        Plots.annotate!(X3[j].+ x[end].+0.1, X2[j]+1.05, Plots.text("-"),dpi=500,size=(1200,800))
        Plots.annotate!(X3[j].+ x[end].+0.25, X2[j]+1.05, Plots.text(ts[j], :left,6),dpi=500,size=(1200,800))
    end

## labelling for the x axis
    for j = 1:length(X1)
        Plots.annotate!(X1[j], 0, Plots.text("|",7),dpi=500,size=(1200,800))
        Plots.annotate!(X1[j], -.4, Plots.text(round(X1[j]),8,:below),dpi=500,size=(1200,800))
    end
    return p
end
