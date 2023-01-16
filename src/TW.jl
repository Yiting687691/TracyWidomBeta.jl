function TW(β;cheb=10^3,method="finite",interp=true,pdf=true,x0=13.0,xN=-10.0,Δx_f=-0.01,Δx_s=-0.002,M_f=Int(floor(-10/(Δx_f))),M_s=Int(floor(-16/(Δx_s))),l=20)
    if method=="finite" && pdf=false && interp=true
        TW=fd(β;cheb,interp=true,pdf=false,x0,xN,Δx_f,M_f)
    elseif method=="finite" && pdf=true && interp=true
        TW=fd(β;cheb,interp=true,pdf=true,x0,xN,Δx_f,M_f)
    elseif method=="finite" && pdf=true && interp=false
        TW=fd(β;cheb,interp=false,pdf=true,x0,xN,Δx_f,M_f)
    elseif method=="finite" && pdf=false && interp=false
        TW=fd(β;cheb,interp=false,pdf=false,x0,xN,Δx_f,M_f)
    elseif method=="bdf3" && pdf=true && interp=true
        TW=sd(β;method=="bdf3",cheb,interp=true,pdf=true,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf3" && pdf=false && interp=true
        TW=sd(β;method=="bdf3",cheb,interp=true,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf3" && pdf=false && interp=false
        TW=sd(β;method=="bdf3",cheb,interp=false,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf3" && pdf=true && interp=false
        TW=sd(β;method=="bdf3",cheb,interp=false,pdf=true,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf4" && pdf=true && interp=true
        TW=sd(β;method=="bdf4",cheb,interp=true,pdf=true,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf4" && pdf=false && interp=true
        TW=sd(β;method=="bdf4",cheb,interp=true,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf4" && pdf=false && interp=false
        TW=sd(β;method=="bdf4",cheb,interp=false,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf4" && pdf=true && interp=false
        TW=sd(β;method=="bdf4",cheb,interp=false,pdf=true,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf5" && pdf=true && interp=true
        TW=sd(β;method=="bdf5",cheb,interp=true,pdf=true,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf5" && pdf=false && interp=true
        TW=sd(β;method=="bdf5",cheb,interp=true,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf5" && pdf=false && interp=false
        TW=sd(β;method=="bdf5",cheb,interp=false,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf5" && pdf=true && interp=false
        TW=sd(β;method=="bdf5",cheb,interp=false,pdf=true,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf6" && pdf=true && interp=true
        TW=sd(β;method=="bdf6",cheb,interp=true,pdf=true,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf6" && pdf=false && interp=true
        TW=sd(β;method=="bdf6",cheb,interp=true,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf6" && pdf=false && interp=false
        TW=sd(β;method=="bdf6",cheb,interp=false,pdf=false,x0,xN,Δx_s,M_s,l)
    elseif method=="bdf6" && pdf=true && interp=false
        TW=sd(β;method=="bdf6",cheb,interp=false,pdf=true,x0,xN,Δx_s,M_s,l)
    else
        TW="Input valid arguments"
    end
    return TW
end
