require(statmod)


# base functions

closedIntervalWeightedSum <- function(f, lower, upper, quad) {

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

halfOpenIntervalWeightedSum <- function(f, finite.limit, quad) {

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

openIntervalWeightedSum <- function(f, quad) {
  sum(f(quad$nodes) * quad$weights)
}


# wrappers

legendreIntegrate <- function(f, lower, upper, subdivisions=16) {
  quad <- statmod::gauss.quad(subdivisions, kind='legendre')
  closedIntervalWeightedSum(f, lower, upper, quad)
}

jacobiIntegrate <- function(f, lower, upper, subdivisions=16, alpha=0, beta=0, divide.weight=TRUE) {
  if (alpha <= -1 || beta <= -1)
    stop("alpha and beta parameters for Jacobi must be > -1")
  quad <- statmod::gauss.quad(subdivisions, kind='jacobi', alpha=alpha, beta=beta)
  if (divide.weight)
    integrand <- function(x) f(x) / ((1-x)^alpha * (1+x)^beta)
  else
    integrand <- f
  closedIntervalWeightedSum(integrand, lower, upper, quad)
}

chebyshev1Integrate <- function(f, lower, upper, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='chebyshev1')
  if (divide.weight)
    integrand <- function(x) f(x) * sqrt(1-x^2)
  else
    integrand <- f
  closedIntervalWeightedSum(integrand, lower, upper, quad)
}

chebyshev2Integrate <- function(f, lower, upper, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='chebyshev2')
  if (divide.weight)
    integrand <- function(x) f(x) / sqrt(1-x^2)
  else
    integrand <- f
  closedIntervalWeightedSum(integrand, lower, upper, quad)
}

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
  ifelse(lower < upper, 1, -1) * halfOpenIntervalWeightedSum(integrand, finite.limit, quad)
}

hermiteIntegrate <- function(f, subdivisions=16, divide.weight=TRUE) {
  quad <- statmod::gauss.quad(subdivisions, kind='hermite')
  if (divide.weight)
    integrand <- function(x) f(x) * exp(x^2)
  else
    integrand <- f
  openIntervalWeightedSum(integrand, quad)
}