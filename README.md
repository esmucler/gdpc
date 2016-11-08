
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/esmucler/gdpc.svg?branch=master)](https://travis-ci.org/esmucler/gdpc)

gdpc
====

This package (currently under development!) provides functions for computing the Generalized Dynamic Principal Components proposed in Peña and Yohai (2016). The number of components can be supplied by the user or chosen automatically so that a given proportion of variance is explained. The number of lags is chosen automatically using one of the following criteria: Leave-one-out cross-validation, an AIC type criterion, a BIC type criterion or a criterion based on a proposal of Bai and Ng (2002).

------------------------------------------------------------------------

### Installation

You can install the **development** version from [GitHub](https://github.com/esmucler/gdpc)

``` r
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

This package is free and open source software, licensed under GPL (&gt;= 2).
