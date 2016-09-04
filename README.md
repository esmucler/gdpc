
<!-- README.md is generated from README.Rmd. Please edit that file -->
### gdpc

This package (currently under development!) provides functions for computing the Generalized Dynamic Principal Components proposed in Peña and Yohai (2016). The number of components can be supplied by the user or chosen automatically so that a given proportion of variance is explained. The number of lags is chosen automatically using one of the following criteria: Leave-one-out cross-validation, an AIC type criterion, a BIC type criterion or a criterion based on a proposal of Bai and Ng (2002).

### Installation

You can install the **development** version from [GitHub](https://github.com/esmucler/gdpc)

``` r
devtools::install_github("esmucler/gdpc")
```

### Usage

``` r
library(gdpc)
#Load and plot Industrial Production Index Data
data(ipi91)
plot(ipi91, plot.type = 'multiple', main = 'Industrial Production Index')
```

![](README-unnamed-chunk-2-1.png)

``` r
#Choose number of lags among 0, \dots, 9 using the Leave One Out criterion
gdpc_ipi <- auto.gdpc(ipi91, k_max = 9, niter_max = 1500, ncores = 4)
#> [1] "Computing component number 1"
#> [1] "Total number of computed components: 1"
gdpc_ipi
#>             Number.of.lags   LOO   MSE Explained.Variance
#> Component 1              9 8.471 7.778              0.948
#Plot the component
plot(gdpc_ipi, which_comp = 1, ylab = '')
```

![](README-unnamed-chunk-2-2.png)

``` r
#Get fitted values and plot
fit_val <- fitted(gdpc_ipi[[1]])
colnames(fit_val) <- colnames(ipi91)
plot(fit_val, main = 'Fitted values')
```

![](README-unnamed-chunk-2-3.png)

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
