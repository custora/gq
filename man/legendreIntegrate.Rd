\name{legendreIntegrate}
\alias{legendreIntegrate}
\title{Computes a Gauss-Legendre quadrature, which has weight function \eqn{w(x) = 1}.
This is used for integration over finite limits.}
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

