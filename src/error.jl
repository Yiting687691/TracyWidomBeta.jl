finite=finite_cdf(2);
F=finite[1];
spec=BDF4_cdf(2);
G=spec[1];
initial_time=8;
delta_x=-0.05;
final_time=-8;
time=initial_time:delta_x:final_time;
T=TracyWidom;
err=zeros(length(time),1);
for j=1:length(time)
    t=time[j];
    err[j]=log10(abs(F(t)-cdf(T,t; beta=2,num_points=300)))
end
err=vec(err);
err2=zeros(length(time),1);
for k=1:length(time)
    t=time[k];
    err2[k]=log10(abs(G(t)-cdf(T,t; beta=2,num_points=300)))
end
err2=vec(err2);
p1=plot(err, label = "Finite Difference Method", lw = 3,aspect_ratio=5,ylims=(-18,-5))
plot!(err2, label = "Spectral Method", lw = 3,aspect_ratio=5,ylims=(-18,-5))
ylabel!("Error")




p2=plot(F', label = "Largest", lw = 3)
plot!(FF', label = "Second Largest", lw = 3)
plot!(FFF', label = "Third Largest", lw = 3)
xlabel!("s")
