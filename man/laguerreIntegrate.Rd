\name{laguerreIntegrate}
\alias{laguerreIntegrate}
\title{Computes a Gauss-Laguerre quadrature}
\usage{
laguerreIntegrate(f, lower, upper, subdivisions = 16, alpha = 0,
  divide.weight = TRUE)
}
\arguments{
  \item{f}{The function to integrate.}

  \item{lower}{The lower limit(s) of integration. May be a
  vector. Either this or \code{upper}, but not both, should
  be infinite.}

  \item{upper}{The upper limit(s) of integration. May be a
  vector. Either this or \code{lower}, but not both should
  be infinite.}

  \item{subdivisions}{The number of terms to use in the
  quadrature sum.}

  \item{alpha}{The \eqn{alpha} used in the weight
  function.}

  \item{divide.weight}{Should f be divided by the weight
  function?}
}
\value{
A vector of computed sums that approximates an integral
over the supplied limits.
}
\description{
Computes a Gauss-Laguerre quadrature, which has weight
function \eqn{w(x) = x^\alpha exp(-x)}. This is used for
integration from a finite limit to positive or negative
infinity, or vice versa.
}
\examples{
laguerreIntegrate(function(x) exp(-x), 0:3, Inf)
laguerreIntegrate(function(x) exp(-x), c(0,1,Inf,Inf), c(Inf,Inf,0,1))
# next example: some NAs and a warning, but the limits that work will produce values
laguerreIntegrate(function(x) exp(-x), c(0,1,Inf,1), c(Inf,2,Inf,Inf))
laguerreIntegrate(function(x) exp(x), -Inf, 0)
}
\references{
http://en.wikipedia.org/wiki/Gaussian_quadrature
}

