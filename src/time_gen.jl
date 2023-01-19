function time_gen(x1,x2,Δx)
    xs=x1:Δx:x2 |> Array
    if abs(xs[end])<abs(x2)
        xs=x1:Δx:(x2+Δx) |> Array
    end
    return xs
end
