\name{components}
\alias{components}
\title{Generic Function for Getting Components From an Object}
\description{
  Generic function for getting components from an object.
}


\usage{
  components(object, which_comp)
}
\arguments{
  \item{object}{An object of class \code{gdpcs}, usually the result of \link{auto.gdpc}.}
  \item{which_comp}{Numeric vector indicating which components to get. Default is 1.}
}


\value{ A matrix whose columns are the desired components.}
\author{
  Daniel Peña, Ezequiel Smucler, Victor Yohai
}



\seealso{
  \code{\link{gdpc}}, \code{\link{auto.gdpc}}, \code{\link{plot.gdpc}}
}


\keyword{ts}