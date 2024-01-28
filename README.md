Compute the cumulative distribution function (CDF) and probability density function (PDF) of the Tracy-Widom distribution for arbitrary β>0.

This respository contains the code for: https://doi.org/10.3842/SIGMA.2024.005

A basic example usage follows.
```julia
β = 3 #value of β
F = TW(β) #CDF
f = TW(β; pdf = true) #PDF
s = -2 #evaluation point
F(s)
f(s)
```
The provided code employs the default finite-difference discretization paired with the trapezoidal method. By setting the 'pdf' argument to true, it can return the PDF of the Tracy-Widom distribution for β=2 and then evaluate the function at x=0.5. Both the CDF and PDF can be seamlessly plotted for x values ranging from -10 to 13/√β. Additionally, the PDF can be acquired using the subsequent approach.
```julia
F_pdf=F_cdf'
```
For enhanced accuracy, you can configure the arguments as shown below. Here, I've opted to apply BDF5 to the spectral discretization. 
```julia
F_cdf=TW(2;method="spectral",step="bdf5")
```
Please be aware that for β=4, x should be divided by a factor of 2^(1/6), as many published results adhere to a distinct scaling convention. For instance, to calculate the cdf of the Tracy-Widom distribution for β=4 at x=-2 using the old scaling convention, refer to the provided code below.
```julia
F_cdf=TW(4)
F_cdf(-2/(2^(1/6)))
```
Also, be aware that both the default finite-difference discretization and the more precise spectral discretization are effective for 1≤β≤30. For values of β greater than 30, it's necessary to employ a finer grid and increase the value of M. For 0<β<1, it's essential to expand the domain for x.

By default, all the provided examples utilize Fourier interpolation to produce a continuous function. However, one can also obtain just the discrete values prior to interpolation. For instance, the code below demonstrates how to retrieve discrete values of the cdf using the finite-difference discretization with the trapezoidal method, along with the associated x values.
```julia
F_cdf=TW(2;interp=false)
x=F[1]
F_cdf=F[2]
```

See [Index](https://github.com/Yiting687691/TracyWidomBeta.jl/blob/main/notebook/Index.ipynb) for a brief summary of how this algorithm is created and related plots.

## References

T. Trogdon & Y. Zhang (2023), Computing the Tracy-Widom distribution for arbitrary β>0, [arXiv:2304.04951](https://arxiv.org/abs/2304.04951)





[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/main/graph/badge.svg?token=Q9ZOX49RPV)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7879558.svg)](https://doi.org/10.5281/zenodo.7879558)
