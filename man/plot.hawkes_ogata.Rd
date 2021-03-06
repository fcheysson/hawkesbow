% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hawkes_ogata.R
\name{plot.hawkes_ogata}
\alias{plot.hawkes_ogata}
\title{Plot of a simulated Hawkes process}
\usage{
\method{plot}{hawkes_ogata}(x, precision = 1000, ...)
}
\arguments{
\item{x}{A simulated Hawkes process from \code{hawkes_ogata}.}

\item{precision}{(default = 1e3) Number of points to plot.}

\item{...}{Only there to fit the declaration of S3 method \code{plot}.}
}
\value{
None
}
\description{
Plots a Hawkes process simulated by the function \code{hawkes_ogata},
highlighting the steps used in Ogata's thinning algorithm.
}
\examples{
# Simulate an exponential Hawkes process with baseline intensity 1 and
# excitation function 1*exp(-2t)
x <- hawkes_ogata(10, 1, 1, 2)
plot(x)
}
