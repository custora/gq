
require(statmod)

# base helper functions, not exported

finiteIntervalWeightedSum <- function(f, lower, upper, quad) {

  lower <- as.numeric(lower)
  upper <- as.numeric(upper)

  if (any(is.na(lower)) || any(is.na(upper))) {
    warning("NAs found in limits")
  }

  if (any(is.infinite(lower)) || any(is.infinite(upper))) {
    warning("infinite limits passed to a closed interval sum")
    lower[is.infinite(lower)] <- NA
    upper[is.infinite(upper)] <- NA
  }

  mids <- (upper+lower)/2
  halfwidths <- (upper-lower)/2
  
  # for loop was actually faster than mapply for typical subdivisions and long 
  # upper and lower limits

  result <- rep(0, max(length(lower), length(upper)))
  for (i in 1:length(quad$nodes))
    result <- result + f(mids + halfwidths * quad$nodes[i]) * quad$weights[i] * halfwidths

  result
}

halfInfiniteIntervalWeightedSum <- function(f, finite.limit, quad) {

  if (any(is.na(finite.limit))) {
    warning("NAs found in limits")
  }

  nonconforming <- is.infinite(finite.limit)
  if (any(nonconforming)) {
    warning("infinite finite.limit limit passed to a half open interval sum")
    finite.limit[nonconforming] <- NA
  }

  finite.limit <- as.numeric(finite.limit)
  result <- rep(0, length(finite.limit))
  for (i in 1:length(quad$nodes))
    result <- result + f(quad$nodes[i]) * quad$weights[i]

  result
}

infiniteIntervalWeightedSum <- function(f, quad) {
  sum(f(quad$nodes) * quad$weights)
}


#' Computes a Gauss-Legendre quadrature, which has weight function \eqn{w(x) = 1}. 
#' This is used for integration over finite limits. 
#'
#' @param f The function to integrate.
#' @param lower The lower limit(s) of integration. May be a vector. Should be 
#'   finite. 
#' @param upper The upper limit(s) of integration. May be a vector. Should be
#'   finite. 
#' @param subdivisions The number of terms to use in the quadrature sum. 
#' @param divide.weight Should f be divided by the weight function? No effect
#'   for Legendre quadrature, since the weight function is 1, but provided
#'   just for consistency with the other integration functions. 
#' @return A vector of computed sums that approximates an integral over the 
#'   supplied limits. 
#' @export
legendreIntegrate <- function(f, lower, upper, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='legendre')
  finiteIntervalWeightedSum(f, lower, upper, quad)
}

#' Computes a Gauss-Jacobi quadrature, which has weight function 
#' \eqn{w(x) = (1-x)^\alpha (1+x)^\beta}. This is used for integration over 
#' finite limits. 
#'
#' @param f The function to integrate.
#' @param lower The lower limit(s) of integration. May be a vector. Should be
#'   finite. 
#' @param upper The upper limit(s) of integration. May be a vector. Should be
#'   finite. 
#' @param subdivisions The number of terms to use in the quadrature sum. 
#' @param alpha The \eqn{alpha} 
#' @param divide.weight Should f be divided by the weight function? 
#' @return A vector of computed sums that approximates an integral over the 
#'   supplied limits. 
#' @export
jacobiIntegrate <- function(f, lower, upper, subdivisions=16, alpha=0, beta=0, divide.weight=TRUE) {
  if (alpha <= -1 || beta <= -1)
    stop("alpha and beta parameters for Jacobi must be > -1")
  quad <- statmod::gauss.quad(subdivisions, kind='jacobi', alpha=alpha, beta=beta)
  if (divide.weight)
    integrand <- function(x) f(x) / ((1-x)^alpha * (1+x)^beta)
  else
    integrand <- f
  finiteIntervalWeightedSum(integrand, lower, upper, quad)
}

#' Computes a Gauss-Chebyshev quadrature of type 1, which has weight function
#' \eqn{w(x) = 1/\sqrt{1-x^2}}. This is used for integration over finite 
#' limits.
#'
#' @param f The function to integrate.
#' @param lower The lower limit(s) of integration. May be a vector. Should be
#'   finite. 
#' @param upper The upper limit(s) of integration. May be a vector. Should be
#'   finite. 
#' @param subdivisions The number of terms to use in the quadrature sum. 
#' @param divide.weight Should f be divided by the weight function? 
#' @return A vector of computed sums that approximates an integral over the 
#'   supplied limits. 
#' @export
chebyshev1Integrate <- function(f, lower, upper, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='chebyshev1')
  if (divide.weight)
    integrand <- function(x) f(x) * sqrt(1-x^2)
  else
    integrand <- f
  finiteIntervalWeightedSum(integrand, lower, upper, quad)
}

#' Computes a Gauss-Chebyshev quadrature of type 2, which has weight function
#' \eqn{w(x) = \sqrt{1-x^2}}. This is used for integration over finite limits. 
#'
#' @param f The function to integrate.
#' @param lower The lower limit(s) of integration. May be a vector. Should be
#'   finite. 
#' @param upper The upper limit(s) of integration. May be a vector. Should be
#'   finite. 
#' @param subdivisions The number of terms to use in the quadrature sum. 
#' @param divide.weight Should f be divided by the weight function? 
#' @return A vector of computed sums that approximates an integral over the 
#'   supplied limits. 
#' @export
chebyshev2Integrate <- function(f, lower, upper, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='chebyshev2')
  if (divide.weight)
    integrand <- function(x) f(x) / sqrt(1-x^2)
  else
    integrand <- f
  finiteIntervalWeightedSum(integrand, lower, upper, quad)
}

#' Computes a Gauss-Laguerre quadrature, which has weight function 
#' \eqn{w(x) = x^\alpha exp(-x)}. This is used for integration from a finite 
#' limit to positive or negative infinity, or vice versa.
#'
#' @param f The function to integrate.
#' @param lower The lower limit(s) of integration. May be a vector. Either this
#'   or \code{upper}, but not both, should be infinite. 
#' @param upper The upper limit(s) of integration. May be a vector. Either this
#'   or \code{lower}, but not both should be infinite. 
#' @param subdivisions The number of terms to use in the quadrature sum. 
#' @param divide.weight Should f be divided by the weight function? 
#' @return A vector of computed sums that approximates an integral over the 
#'   supplied limits. 
#' @export
laguerreIntegrate <- function(f, lower, upper, subdivisions=16, alpha=0, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='laguerre', alpha=alpha)
  lower.finite   <- is.finite(lower)
  lower.infinite <- is.infinite(lower)
  upper.finite   <- is.finite(upper)
  upper.infinite <- is.infinite(upper)
  finite.limit   <- ifelse(lower.finite & upper.infinite, lower,
                      ifelse(lower.infinite & upper.finite, upper, NA))
  sign.inf <- sign(ifelse(lower.infinite, lower, upper))
  g <- function(x) f(sign.inf * (x + finite.limit))
  if (divide.weight)
    integrand <- function(x) g(x) / (x^alpha * exp(-x))
  else
    integrand <- g
  ifelse(lower < upper, 1, -1) * halfInfiniteIntervalWeightedSum(integrand, finite.limit, quad)
}

#' Computes a Gauss-Hermite quadrature, which has weight function 
#' \eqn{w(x) = exp(-x^2)}. This is used for integration from negative to
#' positive infinity. 
#'
#' @param f The function to integrate.
#' @param subdivisions The number of terms to use in the quadrature sum. 
#' @param divide.weight Should f be divided by the weight function? 
#' @return A vector of computed sums that approximates an integral over 
#'   negative to positive infinity. 
#' @export
hermiteIntegrate <- function(f, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='hermite')
  if (divide.weight)
    integrand <- function(x) f(x) * exp(x^2)
  else
    integrand <- f
  infiniteIntervalWeightedSum(integrand, quad)
}