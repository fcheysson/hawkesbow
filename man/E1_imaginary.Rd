% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{E1_imaginary}
\alias{E1_imaginary}
\title{Exponential integral of imaginary argument}
\usage{
E1_imaginary(x)
}
\arguments{
\item{x}{A non-negative number}
}
\value{
The exponential integral of argument ix
}
\description{
Calculates the value of
\deqn{E_1(ix) = \int_1^\infty \frac{e^{-ixt}}{t} \mathrm{d}t}
using its relation to the trigonometric integrals
(cf. \url{https://en.wikipedia.org/wiki/Exponential_integral#Exponential_integral_of_imaginary_argument}):
\deqn{E_1(ix) = i \left[ -\frac{1}{2} \pi + Si(x) \right] - Ci(x)}
and their Pad\'e approximants
(cf. \url{https://en.wikipedia.org/wiki/Trigonometric_integral#Efficient_evaluation})
}
\examples{
E1_imaginary(1.0)
}
