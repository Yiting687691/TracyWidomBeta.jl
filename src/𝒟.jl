function 𝒟(L,m,j)
  diffvec = (L,m,j) -> ((-floor(m/2):1:floor((m-1)/2))*(2*1im*pi/L)).^j
  spdiagm(diffvec(L,m,j))
end
