#This function gives the waterfall plot of the density of the Tracy-Widom distribution function for the values of beta defined by B. Ex: B=1:1:10.
function water_TW(B)
    x=-6:0.01:5
    dd = zeros(length(B),length(x))
    for i=1:length(B)
        beta=B[i]
        F=TW(beta)
        f=F'
        for j=1:length(x)
            dd[i,j]=f(x[j])
        end
    end

    δ= 0.8
    m = length(B)
    slope = 4
    dx = 1
    dt = 1
    sc = 8
    up = 10
    C = 0.1

    X1 = x[1]:dx:x[end]
    X2 = B.-B[1]
    tick_shift = 0.05
    axis_shift = 0.05
    X2 = X2/X2[end]*(δ*(m-1)).+ tick_shift
    X3 = X2/slope
    x_range = [x[1],x[end]+ δ/slope*m]

    p=Plots.plot()
    for j = m:-1:1
        Plots.plot!(x .+ δ*(j-1)/slope, sc*dd[j,:] .+ δ*(j-1) .+ C, grid = false, legend = false, lw = 2, color = :black, fillrange = δ*(j-1), fillcolor = :white, seriesalpha = 2, alpha = 0.9,dpi=1000)
    end

    ## draws the x axis
    Plots.plot!(X1, 0*X1, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]),dpi=1000 )

    ## draws the t axis
    Plots.plot!(X3 .+ x[end].+0.1, X2.+axis_shift, color = :black, yticks = false, xticks = false, xaxis = (false,x_range), yaxis = (false,[0,δ*(m+5) + up]),dpi=1000 )

    ## labelling for the t axis
    for j = 1:length(X3)
        Plots.annotate!(X3[j].+ x[end].+0.1, X2[j]+0.1, Plots.text("-"),dpi=1000)
        Plots.annotate!(X3[j].+ x[end].+0.2, X2[j]+0.1, Plots.text(B[j], :left,6),dpi=1000)
    end

    ## labelling for the x axis

    for j = 1:length(X1)
        Plots.annotate!(X1[j], 0, Plots.text("|",7),dpi=1000)
        Plots.annotate!(X1[j], -.4, Plots.text(round(X1[j]),8,:below),dpi=1000)
    end
    return p
end
