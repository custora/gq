\name{laguerreIntegrate}
\alias{laguerreIntegrate}
\title{Computes a Gauss-Laguerre quadrature, which has weight function
\eqn{w(x) = x^\alpha exp(-x)}. This is used for integration from a finite
limit to positive or negative infinity, or vice versa.}
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

