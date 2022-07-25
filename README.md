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
This code uses the default finite-difference method. It first computes the Tracy-Widom distribution function for β=2, and then evaluates it at x=0.5. The probability density function can also be obtained and be evaluated at any point. Both distributions can be plotted directly. F is defined on the interval [-10,13]. If higher accuracy is preferred, specify the method as the following.
```
F=TW(2,method="BDF4")
f=F'
```
Note that the default method only works well for β⩽450, and the accurate one works well for β⩽10. For larger values of β, one needs to use more Chebyshev points (Both algorithms use 1000 Chebyshev points).


[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![Codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
