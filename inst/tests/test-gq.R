
test_that("Legendre works precisely with polynomials", {

  # basic
  expect_that(legendreIntegrate(function(x) x^2, -1, 1), equals(2/3))
  expect_that(legendreIntegrate(function(x) 5*x^3 - 3*x^2 + 2, -5, 3), equals(-816))

  # reversed limits
  expect_that(legendreIntegrate(function(x) 5*x^3 - 3*x^2 + 2, 3, -5), equals(816))

  # vectorization
  expect_that(legendreIntegrate(function(x) x^2, 0, c(1,2,3)), equals(c(1,2,3)^3/3))
  expect_that(legendreIntegrate(function(x) x^4, c(-4,-2,0), c(1,2,3)), equals((c(1,2,3)^5 - c(-4,-2,0)^5)/5))

  # matrices even if you like
  expect_that(legendreIntegrate(function(x) x^2, matrix(1:9, 3, 3), t(matrix(1:9, 3, 3))), 
              equals((t(matrix(1:9, 3, 3))^3 - matrix(1:9, 3, 3)^3)/3))

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

})

test_that("Laguerre works", {

})

test_that("Hermite works", {

  expect_that(hermiteIntegrate(function(x) exp(-x^2)), equals(sqrt(pi)))
  expect_that(hermiteIntegrate(function(x) exp(-x^2/2)), equals(sqrt(2*pi), tolerance=1e-6))

})

test_that("Jacobi and Chebyshev match in special cases", {

})


test_that("errors on unrecoverably bad inputs", {
  
  expect_that(legendreIntegrate(function(x) x, "a", 0:2), throws_error())

  # infinite legendre or finite other
  expect_that(legendreIntegrate(function(x) x, 0, Inf), throws_error())
  expect_that(jacobiIntegrate(function(x) x, 0, Inf), throws_error())
  expect_that(chebyshev1Integrate(function(x) x, 0, Inf), throws_error())
  expect_that(chebyshev2Integrate(function(x) x, 0, Inf), throws_error())
  expect_that(laguerreIntegrate(function(x) x, -1, 1), throws_error())

})

test_that("warnings on some bad inputs", {
  
  # inputs of mismatched size
  expect_that(legendreIntegrate(function(x) x, 1:3, 1:4), gives_warning())

  # NAs

})

test_that("better matches with ")