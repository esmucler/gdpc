
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/esmucler/gdpc.svg?branch=master)](https://travis-ci.org/esmucler/gdpc)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/gdpc)](https://cran.r-project.org/package=gdpc)
[![Downloads](http://cranlogs.r-pkg.org/badges/gdpc)](https://cran.r-project.org/package=gdpc)
![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/gdpc)

# gdpc

This package provides functions for computing the Generalized Dynamic
Principal Components proposed in [Peña and Yohai
(2016)](http://dx.doi.org/10.1080/01621459.2015.1072542). It also
provides a function to automatically choose the number of components and lags, using
the criteria described in [Peña, Smucler and Yohai
(2020)](http://dx.doi.org/10.18637/jss.v092.c02).

-----

### Installation

You can install the **stable** version from [R
CRAN](https://cran.r-project.org/package=gdpc).

``` r
install.packages('gdpc', dependencies = TRUE)
```

You can install the **development** version from
[GitHub](https://github.com/esmucler/gdpc)

``` r
library(devtools)
devtools::install_github("esmucler/gdpc")
```

### Usage

``` r
library(gdpc)
# An example using an artificial data set
T <- 200 
m <- 500
set.seed(1234)
f <- rnorm(T + 1)
x <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
}
# Compute first DPC using one lag
fit <- gdpc(x, k = 1)
fit
# Plot loadings
par(mfrow = c(1, 2))
plot(fit, which = 'Loadings', which_load = 0, xlab = '', ylab = '') 
plot(fit, which = 'Loadings', which_load = 1, xlab = '', ylab = '') 
```

``` r
# An example using a real data set
# Load and plot Industrial Production Index Data
data(ipi91)
par(mfrow = c(1, 1))
plot(ipi91, plot.type = 'multiple', main = 'Industrial Production Index')
# Choose number of lags among 0,..., 9 using the (default) Leave One Out criterion
# This might take a minute.
gdpc_ipi <- auto.gdpc(ipi91, k_max = 9, niter_max = 1500, ncores = 4)
gdpc_ipi
# Plot the component
plot(gdpc_ipi, which_comp = 1, ylab = '')
# Get reconstruction of the time series and plot
recons <- fitted(gdpc_ipi, 1)
colnames(recons) <- colnames(ipi91)
plot(recons, main = 'Fitted values')
```

### License

This package is free and open source software, licensed under GPL (\>=
2).
