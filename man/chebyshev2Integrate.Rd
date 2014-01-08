\name{chebyshev2Integrate}
\alias{chebyshev2Integrate}
\title{Computes a Gauss-Chebyshev quadrature of type 2, which has weight function
\eqn{w(x) = \sqrt{1-x^2}}. This is used for integration over finite limits.}
\usage{
chebyshev2Integrate(f, lower, upper, subdivisions = 16,
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

  \item{divide.weight}{Should f be divided by the weight
  function?}
}
\value{
A vector of computed sums that approximates an integral
over the supplied limits.
}
\description{
Computes a Gauss-Chebyshev quadrature of type 2, which has
weight function \eqn{w(x) = \sqrt{1-x^2}}. This is used for
integration over finite limits.
}

