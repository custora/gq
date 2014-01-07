
# gq

[Gaussian quadrature](http://en.wikipedia.org/wiki/Gaussian_quadrature)-based integration functions in R. Depends on the [statmod](http://cran.r-project.org/web/packages/statmod/index.html) package (in fact statmod's `gauss.quad` function computes the quadrature points and weights). The integration limits can be vectors. 

There is no numerical optimization to find the exact answer or to get the answer right within a certain error tolerance - you simply specify the number of points you want to use in the approximation. So this is not the best approach when you need high accuracy. This is more for when you want fast computation over a large number of different limits. 

As an example of where you'd want that: we use this to evaluate integrals in likelihood functions. We may have millions of observations, each with a different limit, but in practice it may not adversely affect likelihood maximization too much if the integral is not exactly calculated. 

## Installation

## Usage

## Tests
