# TracyWidomBeta
Compute both the cdf and pdf of the Tracy-Widom distribution for arbitrary β>0.

An example usage follows.
```
F=TW(2,1000)
F_cdf=F[1]
F_cdf(0.5)
F_pdf=F[2]
F_pdf(0.5)
plot(F_cdf)
plot(F_pdf)
```
The above code uses the default finite-difference discretization. It returns both cdf and pdf of the Tracy-Widom distribution function for β=2 using 1000 Fourier modes, and then evaluates them at x=0.5. The pdf can also be obtained by the following way.
```
F_pdf=F_cdf'
```
Both cdf and pdf can be plotted directly for x between -10 and 13. If higher accuracy is preferred, specify the method as the following.
```
F=TW(2,1000;method="spectral")
```
Note that for β=4, one need to account for an extra 2^(1/6) factor when evaluating at some x value. For example, if one wants to compute the cdf of the Tracy-Widom distribution for β=4 at x=-2 using the finite-difference discretization, use the following code.
```
F=TW(4,1000)
F_cdf=F[1]
F_cdf(-2/(2^(1/6)))
```
Also note that the default finite-difference discretization only works well for 1⩽β⩽400, and the more accurate spectral discretization works well for 1⩽β⩽10. For larger values of β, one needs to use finer grid and more Chebyshev points (both methods use 1000 Chebyshev points by default). For β<1, one needs to enlarge the domain for x.


[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![Codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
