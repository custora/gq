\name{jacobiIntegrate}
\alias{jacobiIntegrate}
\title{Computes a Gauss-Jacobi quadrature}
\usage{
jacobiIntegrate(f, lower, upper, subdivisions = 16, alpha = 0, beta = 0,
  divide.weight = TRUE)
}
\arguments{
  \item{f}{The function to integrate.}

  \item{lower}{The lower limit(s) of integration. May be a
  vector. Should be finite.}

  \item{upper}{The upper limit(s) of integration. May be a
  vector. Should be finite.}

  \item{subdivisions}{The number of terms to use in the
  quadrature sum.}

  \item{alpha}{The \eqn{alpha} used in the weight
  function.}

  \item{beta}{The \eqn{beta} used in the weight function.}

  \item{divide.weight}{Should f be divided by the weight
  function?}
}
\value{
A vector of computed sums that approximates an integral
over the supplied limits.
}
\description{
Computes a Gauss-Jacobi quadrature, which has weight
function \eqn{w(x) = (1-x)^\alpha (1+x)^\beta}. This is
used for integration over finite limits.
}
\examples{
jacobiIntegrate(function(x) x^3 + 3*x, 1:3, 2:4) # effectively Legendre
jacobiIntegrate(function(x) sqrt((1+x)/(1-x)), -1, 1, alpha=0, beta=0)       # bad approximation
jacobiIntegrate(function(x) sqrt((1+x)/(1-x)), -1, 1, alpha=-0.5, beta=0.5)  # better
# last element of below will be NA, but the rest will work with a warning
# next example: some NAs and a warning, but the limits that work will produce values
jacobiIntegrate(function(x) x^3 + 3*x, 1:5, c(2:4, Inf, NA))
}
\references{
http://en.wikipedia.org/wiki/Gaussian_quadrature
}

