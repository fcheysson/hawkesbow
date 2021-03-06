% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hawkes.R
\name{hawkes}
\alias{hawkes}
\title{Simulation of a Hawkes process}
\usage{
hawkes(end, fun, repr, family, M = NULL, ...)
}
\arguments{
\item{end}{A non-negative numeric value - right bound of the interval \eqn{[0,\mathrm{end}]}.}

\item{fun}{A non-negative function or numeric value - intensity (function) of the immigrant process.}

\item{repr}{A non-negative numeric value - mean number of offsprings.}

\item{family}{A character string "name" naming a distribution with corresponding random generation function \code{rname}, or directly the random generation function.}

\item{M}{(default = NULL) A non-negative numeric value - upper bound on \code{fun}(ignored if \code{fun} is a numeric value).}

\item{...}{Additional arguments passed on to the random generation function.}
}
\value{
A S3 object of class Hawkes containing a vector ($p) of simulated values,
and all other objects used for the simulation.
}
\description{
Simulates a Hawkes process using its cluster representation:
\itemize{
\item First, the immigrants are drawn according to an (inhomogeneous) Poisson process with intensity measure \code{fun}.
\item Second, the number of offsprings of an immigrant follows a Poisson distribution with intensity \code{repr}.
\item Third, these offsprings are distributed according to the \code{family} distribution.
\item Then, generate further offsprings according to the last two steps.
}
}
\examples{
# Simulate an exponential Hawkes process with baseline intensity 1,
# reproduction mean 0.5 and exponential fertility function with rate 2.
x <- hawkes(10, fun=1, repr=0.5, family="exp", rate=2)
# Simulate a Hawkes process with baseline intensity function 1 + sin(x),
# reproduction mean 0.5 and custom [0,1]-triangular fertility function.
x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5,
            family=function(n) {1 - sqrt(1 - runif(n))})
}
