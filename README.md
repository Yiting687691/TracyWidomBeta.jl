# TracyWidomBeta
Compute the Tracy-Widom distribution for arbitrary β>0.

An example usage follows.
```
F=TW(2)
F(0.5)
f=F'
f(0.5)
plot(F)
plot(f)
```
This code uses the default finite-difference method. It first computes the cdf of the Tracy-Widom distribution function for β=2, and then evaluates it at x=0.5. The pdf can also be obtained and evaluated at x=0.5. Both distributions can be plotted directly. Note that F is defined only on the interval [-10,13]. This domain is so chosen since F is identically zero outside of it. If higher accuracy is preferred, specify the method as the following.
```
F=TW(2,method="BDF4")
f=F'
```
Note that if β=4, one need to account for an extra 2^(1/6) factor when evaluating at some x value. For example, if one wants to compute the cdf of the Tracy-Widom distribution for β=4 at x=-2, use the following codes.
```
F=TW(4)
F(-2/(2^(1/6)))
```
Also note that the default algorithm only works well for β⩽450, and the accurate one works well for β⩽10. For larger values of β, one needs to use more Chebyshev points (Both algorithms use 1000 Chebyshev points by default).


[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![Codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
