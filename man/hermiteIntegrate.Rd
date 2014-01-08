\name{hermiteIntegrate}
\alias{hermiteIntegrate}
\title{Computes a Gauss-Hermite quadrature}
\usage{
hermiteIntegrate(f, subdivisions = 16, divide.weight = TRUE)
}
\arguments{
  \item{f}{The function to integrate.}

  \item{subdivisions}{The number of terms to use in the
  quadrature sum.}

  \item{divide.weight}{Should f be divided by the weight
  function?}
}
\value{
A vector of computed sums that approximates an integral
over negative to positive infinity.
}
\description{
Computes a Gauss-Hermite quadrature, which has weight
function \eqn{w(x) = exp(-x^2)}. This is used for
integration from negative to positive infinity.
}
\examples{
hermiteIntegrate(function(x) exp(-x^2/2))
}
\references{
http://en.wikipedia.org/wiki/Gaussian_quadrature
}

