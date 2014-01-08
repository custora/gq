\name{jacobiIntegrate}
\alias{jacobiIntegrate}
\title{Computes a Gauss-Jacobi quadrature, which has weight function
\eqn{w(x) = (1-x)^\alpha (1+x)^\beta}. This is used for integration over
finite limits.}
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

  \item{alpha}{The \eqn{alpha}}

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

