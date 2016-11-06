\name{gdpc-package}
\alias{gdpc-package}
\docType{package}

\title{Generalized Dynamic Principal Components}

\description{
Computes the Generalized Dynamic Principal Components proposed in Peña and Yohai (2016).
}
\details{
\tabular{ll}{
Package: \tab gdpc\cr
Type: \tab Package\cr
Version: \tab 0.99.3\cr
Date: \tab 2016-11-5\cr
License: \tab GPL (>= 2)\cr
Imports: \tab xts, zoo, Rcpp, parallel, doParallel, foreach\cr
LinkingTo: \tab Rcpp, RcppArmadillo\cr
NeedsCompilation: \tab yes\cr
}

Index:
\preformatted{
auto.gdpc   Computes Generalized Dynamic Principal Components, number of 
            components and lags can be chosen automatically.
}
\preformatted{
components.gdpcs    Get Generalized Dynamic Principal Components from a gdpcs object.
}
\preformatted{
fitted.gdpcs    Get reconstructed time series from a gdpcs object.
}
\preformatted{
gdpc    Computes a single Generalized Dynamic Principal Component with
        a given number of lags.
}
\preformatted{
ipi91   Six series corresponding to the Industrial Production
        Index (IPI) of France, Germany, Italy, United Kingdom, USA and 
        Japan. Monthly data from January 1991 to December 2012.
}
\preformatted{
plot.gdpc   Plots a gdpc object.
}
\preformatted{
plot.gdpcs   Plots a gdpcs object.
}
\preformatted{
pricesSP50   Fifty series corresponding to the stock prices of the first 50
             components of the Standard&Poor's 500 index. Five hundred daily
             observations starting 1/1/2010.
}
}
\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai

Maintainer: Ezequiel Smucler <ezequiels.90@gmail.com>
}

\references{
Peña D. and Yohai V.J. (2016). “Generalized Dynamic Principal Components.” Journal of the
American Statistical Association, 111(515), 1121–1131.
}


\examples{
data(ipi91)
\dontrun{
  #Compute GDPC, number of components and number of lags is chosen automatically.
  #This might take a bit.
  ipi_autogdpc <- auto.gdpc(ipi91)
  }
}
\keyword{ts}
\keyword{package}