# TracyWidomBeta
Compute both the cdf and pdf of the Tracy-Widom distribution for arbitrary β>0.

See (https://arxiv.org/abs/2304.04951)[arxiv.org/abs/2304.04951]

An example usage follows.
```
F_cdf=TW(2)
F_cdf(0.5)
F_pdf=TW(2;pdf=true)
F_pdf(0.5)
plot(F_cdf)
plot(F_pdf)
```
The above code uses the default finite-difference discretization with trapezoidal method. It returns either cdf or pdf of the Tracy-Widom distribution function for β=2, and then evaluates the function at x=0.5. Both cdf and pdf can be plotted directly for x between -10 and 13. The pdf can also be obtained by the following way.
```
F_pdf=F_cdf'
```
If higher accuracy is preferred, specify the arguments as the following, where I choose to apply BDF5 on spectral discretization.
```
F=TW(2;method="spectral",step="bdf5")
```
Note that for β=4, one need to account for an extra 2^(1/6) factor when evaluating at the desired value of x. For example, if one wants to compute the cdf of the Tracy-Widom distribution for β=4 at x=-2 using the finite-difference discretization with trapezoidal method, use the following code.
```
F_cdf=TW(4)
F_cdf(-2/(2^(1/6)))
```
Also note that the default finite-difference discretization only works well for 1⩽β⩽400, and the more accurate spectral discretization works well for 1⩽β⩽10. For larger values of β, one needs to use finer grid and more Chebyshev points (both discretizations use 1000 Chebyshev points by default). For β<1, one needs to enlarge the domain for x.

By default, all of the above examples use Fourier interpolation to generate a continuous function. One can also get only the discrete values before interpolation. For example, the following code shows how to get the discrete values of the cdf using the finite-difference discretization with trapezoidal method and the corresponding values of x.
```
F=TW(2;interp=false)
x=F[1]
F_cdf=F[2]
```

See [Index](https://github.com/Yiting687691/TracyWidomBeta.jl/blob/main/notebook/Index.ipynb) for a brief summary of how this algorithm is created and related plots.


[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![Codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
