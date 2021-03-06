% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hawkes.R
\name{inhpois}
\alias{inhpois}
\title{Simulation of an inhomogeneous Poisson process by thinning}
\usage{
inhpois(end, fun, M = NULL)
}
\arguments{
\item{end}{A non-negative numeric value - right bound of the interval \eqn{[0,\mathrm{end}]}.}

\item{fun}{A non-negative function or numeric value - intensity (function) of the Poisson process.}

\item{M}{(default = NULL) A non-negative numeric value - upper bound on \code{fun} (ignored if \code{fun} is a numeric value).}
}
\value{
A S3 object of class \code{inhpois} containing a vector ($p) of simulated values,
and all other objects used for the simulation.
}
\description{
Simulates an inhomogeneous Poisson process via Ogata's modified thinning algorithm on \eqn{[0,\mathrm{end}]}.
An homogeneous Poisson process with intensity \code{M} is first generated on \eqn{[0,\mathrm{end}]}, then thinned using the specified intensity function \code{fun}.
}
\examples{
# Simulate an inhomogeneous Poisson process with function intensity 1 + sin(x) (bounded by 2)
x <- inhpois(end=10, fun=function(y) {1 + sin(y)}, M=2)
# Simulate a homogeneous Poisson process with intensity 3
x <- inhpois(end=10, fun=3)
}
