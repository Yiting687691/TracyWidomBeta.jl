function TW(β;cheb=10^3,method="finite",interp=true,pdf=false,x0=13.0,xN=-10.0,Δx_f=-0.01,Δx_s=-0.002,M_f=Int(floor(-10/(Δx_f))),M_s=Int(floor(-16/(Δx_s))),l=20)
    if method=="finite"
        TW=fd(β;cheb,interp,pdf,x0,xN,Δx_f,M_f)
    elseif method=="bdf3" || method=="bdf4" || method=="bdf5" || method=="bdf6"
        TW=sd(β;method,cheb,interp,pdf,x0,xN,Δx_s,M_s,l)
    else
        TW="Input valid arguments"
    end
    return TW
end
