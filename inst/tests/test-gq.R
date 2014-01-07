
test_that("Legendre works precisely with polynomials", {

  # basic
  expect_that(legendreIntegrate(function(x) x^2, -1, 1), equals(2/3))
  expect_that(legendreIntegrate(function(x) 5*x^3 - 3*x^2 + 2, -5, 3), equals(-816))

  # reversed limits
  expect_that(legendreIntegrate(function(x) 5*x^3 - 3*x^2 + 2, 3, -5), equals(816))

  # vectorization
  expect_that(legendreIntegrate(function(x) x^2, 0, c(1,2,3)), equals(c(1,2,3)^3/3))
  expect_that(legendreIntegrate(function(x) x^4, c(-4,-2,0), c(1,2,3)), equals((c(1,2,3)^5 - c(-4,-2,0)^5)/5))

})

test_that("Chebyshev1 works", {

  # basic
  expect_that(chebyshev1Integrate(function(x) 1/sqrt(1-x^2), -1, 1), equals(pi))
  expect_that(chebyshev1Integrate(function(x) 1, -1, 1, divide.weight=FALSE), equals(pi))

})

test_that("Chebyshev2 works", {

  # basic
  expect_that(chebyshev2Integrate(function(x) sqrt(1-x^2), -1, 1), equals(0.5 * pi))

})

test_that("Jacobi works", {

  # basic
  expect_that(jacobiIntegrate(function(x) x^2, -1, 1, alpha=0, beta=0), equals(2/3))

  # alpha and beta > -1
  expect_that(jacobiIntegrate(function(x) x, -1, 1, alpha=-2, beta=0), throws_error())
  expect_that(jacobiIntegrate(function(x) x, -1, 1, alpha=0, beta=-2), throws_error())

})

test_that("Laguerre works", {

  # basic
  expect_that(laguerreIntegrate(function(x) exp(-x), 0, Inf), equals(1))
  expect_that(laguerreIntegrate(function(x) exp(-x), Inf, 0), equals(-1))
  expect_that(laguerreIntegrate(function(x) exp(x), -Inf, 0), equals(1))
  expect_that(laguerreIntegrate(function(x) exp(x), 0, -Inf), equals(-1))

  expect_that(laguerreIntegrate(function(x) exp(-x), c(0,1,2), Inf), equals(exp(c(0,-1,-2))))

  expect_that(suppressWarnings(laguerreIntegrate(function(x) exp(-x), c(NA,1,2), c(1,NA,Inf))), equals(c(NA, NA, exp(-2)), tolerance=1e-6))

})

test_that("Hermite works", {

  expect_that(hermiteIntegrate(function(x) exp(-x^2)), equals(sqrt(pi)))
  expect_that(hermiteIntegrate(function(x) exp(-x^2/2)), equals(sqrt(2*pi), tolerance=1e-6))

})

test_that("warnings on some bad inputs", {
  
  # vector inputs of mismatched size
  expect_that(legendreIntegrate(function(x) x, 1:3, 1:4), gives_warning())

  # infinite legendre or finite other
  expect_that(legendreIntegrate(function(x) x, 0, Inf), gives_warning())
  expect_that(jacobiIntegrate(function(x) x, 0, Inf, 1, 1), gives_warning())
  expect_that(chebyshev1Integrate(function(x) x, 0, Inf), gives_warning())
  expect_that(chebyshev2Integrate(function(x) x, 0, Inf), gives_warning())
  expect_that(laguerreIntegrate(function(x) x, -1, 1), gives_warning())

  # NAs are problems and generate warnings, other entries work
  expect_that(legendreIntegrate(function(x) x, NA, 0:2), gives_warning())
  expect_that(suppressWarnings(legendreIntegrate(function(x) x, c(0,NA), c(1,2))), equals(c(0.5, NA)))
  expect_that(legendreIntegrate(function(x) x, c(0,NA), c(1,2)), gives_warning())
  expect_that(legendreIntegrate(function(x) x, "a", 0:2), gives_warning())

})
