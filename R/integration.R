require(statmod)


# base functions, don't call these directly

closedIntervalIntegrate <- function(f, lower, upper, quad, subdivisions=16) {

  if (!is.finite(lower) || !is.finite(upper))
    stop("non-finite limits passed to a closed-interval integration")

  mids <- (upper+lower)/2
  halfwidths <- (upper-lower)/2
  
  # for loop was actually faster than mapply for typical n and long upper and
  # lower vectors

  result <- rep(0, max(length(lower), length(upper)))
  for (i in 1:subdivisions)
    result <- result + f(mids + halfwidths * quad$nodes[i]) * quad$weights[i] * halfwidths

  result
}

halfOpenIntervalIntegrate <- function(f, lower, quad, subdivisions=16) {

  result <- rep(0, length(lower))
  f(quad$nodes[i] + lower)

}

openIntervalIntegrate <- function(f, quad, subdivisions=16) {
  sum(f(quad$nodes) * quad$weights)
}


# wrappers, call these

legendreIntegrate <- function(f, lower, upper, subdivisions=16) {
  quad <- statmod::gauss.quad(subdivisions, kind='legendre')
  closedIntervalIntegrate(f, lower, upper, quad, subdivisions)
}

jacobiIntegrate <- function(f, lower, upper, subdivisions=16, alpha=0, beta=0, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='jacobi')
  if (divide.weight)
    integrand <- function(x) f(x) / ((1-x)^alpha * (1+x)^beta)
  else
    integrand <- f
  closedIntervalIntegrate(integrand, lower, upper, quad, subdivisions, alpha, beta)
}

chebyshev1Integrate <- function(f, lower, upper, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='chebyshev1')
  if (divide.weight)
    integrand <- function(x) f(x) * sqrt(1-x^2)
  else
    integrand <- f
  closedIntervalIntegrate(integrand, lower, upper, quad, subdivisions)
}

chebyshev2Integrate <- function(f, lower, upper, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='chebyshev2')
  if (divide.weight)
    integrand <- function(x) f(x) / sqrt(1-x^2)
  else
    integrand <- f
  closedIntervalIntegrate(integrand, lower, upper, quad, subdivisions)
}

hermiteIntegrate <- function(f, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='hermite')
  if (divide.weight)
    integrand <- function(x) f(x) * exp(x^2)
  else
    integrand <- f
  openIntervalIntegrate(integrand, quad, subdivisions)
}