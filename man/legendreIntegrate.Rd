\name{legendreIntegrate}
\alias{legendreIntegrate}
\title{Computes a Gauss-Legendre quadrature}
\usage{
legendreIntegrate(f, lower, upper, subdivisions = 16, divide.weight = TRUE)
}
\arguments{
  \item{f}{The function to integrate.}

  \item{lower}{The lower limit(s) of integration. May be a
  vector. Should be finite.}

  \item{upper}{The upper limit(s) of integration. May be a
  vector. Should be finite.}

  \item{subdivisions}{The number of terms to use in the
  quadrature sum.}

  \item{divide.weight}{Should f be divided by the weight
  function? No effect for Legendre quadrature, since the
  weight function is 1, but provided just for consistency
  with the other integration functions.}
}
\value{
A vector of computed sums that approximates an integral
over the supplied limits.
}
\description{
Computes a Gauss-Legendre quadrature, which has weight
function \eqn{w(x) = 1}. This is used for integration over
finite limits.
}
\examples{
legendreIntegrate(function(x) x^3 + 3*x, 1:3, 2:4)
legendreIntegrate(function(x) x^3 + 3*x, 2:4, 1:3)
legendreIntegrate(function(x) x^3 + 3*x, 1, 2:4)
# next example: some NAs and a warning, but the limits that work will produce values
legendreIntegrate(function(x) x^3 + 3*x, 1:5, c(2:4, Inf, NA))
}
\references{
http://en.wikipedia.org/wiki/Gaussian_quadrature
}

