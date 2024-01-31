Compute the cumulative distribution function (CDF) and probability density function (PDF) of the Tracy-Widom distribution for arbitrary β>0.

This respository contains the code for: https://doi.org/10.3842/SIGMA.2024.005

A basic example usage follows.
```julia
β = 3 #value of β
F = TW(β) #CDF
f = TW(β; pdf = true) #PDF
s = -2 #point of evaluation
F(s)
f(s)
```
Optional arguments include the method of discretization and the method of time integration. Example usage of these optional arguments follows.
```julia
β = 3 #value of β
F = TW(β; method = "spectral", step = "bdf5") #CDF
f = TW(β; method = "spectral", step = "bdf5", pdf = true) #PDF
s = -2 #point of evaluation
F(s)
f(s)
```
Additionally, discrete values obtained prior to interpolation are accessible. Example usage follows.
```julia
β = 3 #value of β
F = TW(β; interp = false)
x = F[1] #set of discrete points
F_cdf = F[2] #CDF evaluated at these discrete points
```
For β = 4, s should be divided by an extra factor of 2^(1/6), as many published results adhere to a distinct scaling convention.


[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/main/graph/badge.svg?token=Q9ZOX49RPV)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7879558.svg)](https://doi.org/10.5281/zenodo.7879558)
