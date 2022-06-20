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
This code uses the finite-difference method by default and plots the probability density function of the Tracy-Widom distribution. F is defined on the interval [-10,13]. If higher accuracy is preferred, specify the method as follows.
```
F=TW(2,method="BDF4")
f=F'
```

[![CI](https://github.com/Yiting687691/TracyWidomBeta.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Yiting687691/TracyWidomBeta.jl/actions)
[![Codecov](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Yiting687691/TracyWidomBeta.jl)
