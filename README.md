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
This code uses the finite-difference method by default. It first computes the Tracy-Widom distribution function for β=2, and then evaluates it at x=0.5. The probability density function can also be obtained and be evaluated at any point. Both distributions can be plotted directly. F is defined on the interval [-10,13]. If higher accuracy is preferred, specify the method as follows.
```
F=TW(2,method="BDF4")
f=F'
```
Note that the default method only works well for β⩽10^3, for larger values of β, one needs to have finer grid spacing (the default grid spacing is set to be Δx = -0.001):
```
F=TW(10^4,delta_x = -0.001/100)
```
This is also the case for more accurate method, which only works well for β⩽50. One needs to specify not only the method, but also finer grid spacing. It's recommended to use the default method if β is large, otherwise, it'll take forever to get something meaningful.


[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![Codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
