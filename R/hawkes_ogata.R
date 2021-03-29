##############################################################################
#                        SIMULATION OF HAWKES PROCESS                        #
##############################################################################

# Hawkes process by thinning (Ogata's modified thinning algorithm)

#' Simulation of a Hawkes process
#'
#' Simulates a Hawkes process via Ogata's modified thinning algorithm on \eqn{[0,\mathrm{end}]}.
#' This is less efficient than function `hawkes`, but can be useful for pedagogical reasons.
#'
#' @param end Right bound on time.
#' @param lambda Baseline intensity.
#' @param alpha Parameter for the amplitude of the spike.
#' @param beta Parameter for the speed of exponential decay.
#' @param lambda0 (Optional) Initial value of the conditional intensity.
#'
#' @return A S3 object of class Hawkes containing a vector ($p) of simulated values,
#' and all other objects used for the simulation.
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1 and
#' # excitation function 1*exp(-2t)
#' x <- hawkes_ogata(10, 1, 1, 2)
#' plot(x)
hawkes_ogata <- function(end, lambda, alpha, beta, lambda0=NULL) {

	# If no lambda0 specified, take lambda
	# NOTE: in the future, find the correct distribution for lambda0
	if (is.null(lambda0)) lambda0 <- lambda

	# lambda0 must be higher than lambda
	# (else Ogata's algorithm does not work)
	else if (lambda0 < lambda)
		stop("For Ogata's algorithm, lambda0 must be higher than lambda.")

	# Object to return
	sim <- list(p=numeric(0), end=end,
			lambda=lambda, alpha=alpha, beta=beta, lambda0=lambda0,
			n=0, A=numeric(0), M=lambda,
			t=numeric(0), e=numeric(0), U=numeric(0), m=0)
	class(sim) <- "hawkes_ogata"

	# We do as if there was a point in t = 0 that we will
	# remove at the end of the algorithm
	M <- lambda0
	A <- (lambda0 - lambda) / alpha - 1
	p <- 0
	e <- NA
	t <- 0
	U <- NA

	# Initialize
	i <- 1
	n <- 1
	M[2] <- lambda0

	# Ogata's modified thinning algorithm
	while (t[n] < end) {
		n <- n + 1
		e[n] <- rexp(1, rate=M[n])
		t[n] <- t[n-1] + e[n]
		U[n] <- runif(1, min=0, max=M[n])

		if (t[n] < end && U[n] < (lambda + alpha * (A[i]+1) * exp(-beta*(t[n]-p[i])))) {
			A[i+1] <- exp(-beta*(t[n]-p[i]))*(1+A[i])
			p[i+1] <- t[n]
			M[n+1] <- lambda + alpha*(A[i+1]+1)
			i <- i + 1
		} else {
			M[n+1] <- lambda + (M[n] - lambda) * exp(-beta * e[n])
		}
	}

	# Remove the first artificial point
	sim$p <- p[-1]
	sim$M <- M[-1]
	sim$A <- A[-1]
	sim$t <- t[-1]
	sim$e <- e[-1]
	sim$U <- U[-1]
	sim$n <- i-1
	sim$m <- n-1

	return( sim )
}

#' Plot of a simulated Hawkes process
#'
#' Plots a Hawkes process simulated by the function `hawkes_ogata`,
#' highlighting the steps used in Ogata's thinning algorithm.
#'
#' @param x A simulated Hawkes process from `hawkes_ogata`.
#' @param precision (default = 1e3) Number of points to plot.
#' @param ... Only there to fit the declaration of S3 method `plot`.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1 and
#' # excitation function 1*exp(-2t)
#' x <- hawkes_ogata(10, 1, 1, 2)
#' plot(x)
plot.hawkes_ogata <- function(x, precision=1e3, ...) {
	# Conditional intensity
	matplot(z <- seq(0, x$end, by=x$end / precision),
		intensity_ogata(x, z),
		type="l", ylim=c(0, max(x$M)),
		xlab=expression(italic(t)), ylab=expression(italic(U)))
	# Upper bound M of Ogata's modified thinning algorithm
	segments(x0=c(0, x$t[-x$m]), x1=c(x$t[-x$m], x$end),
		y0=x$M[-1-x$m],
		lwd=4, col="blue")
	# All points considered and the corresponding value for U
	points(x=x$t[-x$m],
		y=x$U[-x$m],
		pch=ifelse(x$U[-x$m] > sapply(x$t[-x$m], function(t) {intensity_ogata(x, t)}), 3, 1),
		col=ifelse(x$U[-x$m] > sapply(x$t[-x$m], function(t) {intensity_ogata(x, t)}), "firebrick1", "green4"))
	# The resulting points of the Hawkes process
	points(x=x$p, y=rep(0, x$n), col="green4", pch=15)
	# Nice lines showing which point considered resulted in which point of the Hawkes process
	segments(x0=x$p, y0=rep(0, x$n),
		y1=(x$U[-x$m])[which(x$U[-x$m] <= sapply(x$t[-x$m], function(t) {intensity_ogata(x, t)}))],
		col="green4", lty=2)
	# Legends
	legend("topleft", legend=c(expression(italic(lambda(t))), expression(italic(M)), expression(Accepted), expression(Rejected)),
		lty=c(1,1,NA,NA), col=c("black", "blue", "green", "red"), lwd=c(1,4,NA,NA), pch=c(NA, NA, 1, 3))
}

# #' Intensity of a simulated Hawkes process
# #'
# #' Outputs the intensity of a simulated Hawkes process.
# #'
# #' @param hawkes A simulated Hawkes process from `hawkes_ogata`.
# #' @param t A numeric or vector.
# #'
# #' @return The intensity at time t.
# #'
# #' @export
# #'
# #' @examples
# #' # Simulate an exponential Hawkes process with baseline intensity 1 and
# #' # excitation function 1*exp(-2t)
# #' x <- hawkes_ogata(10, 1, 1, 2)
# #' plot(x)
# #' intensity_ogata(x, 0:10)
intensity_ogata <- function(hawkes, t) {
	# If outside of bounds, return error
	if (any(t < 0) || any(t > hawkes$end))
		stop("t is out of bound.")
	# Create vector of past closest points of hawkes
	index <- sapply(t, function(j) {
		Position(function(p) p >= j, hawkes$p, nomatch=hawkes$n+1)-1
	})
	int <- rep(hawkes$lambda, length(t))
	pind <- index > 0
	int[!pind] <- int[!pind] +
		(hawkes$lambda0 - hawkes$lambda) * exp(-hawkes$beta * t[!pind])
	int[pind] <- int[pind] +
		hawkes$alpha * (hawkes$A[index[pind]]+1) *
		exp(- hawkes$beta * (t[pind] - hawkes$p[index[pind]]))
	return( int )
}

# #' Compensator of a simulated Hawkes process
# #'
# #' Outputs the compensator (integrated intensity) of a simulated Hawkes process.
# #'
# #' @param hawkes A simulated Hawkes process from `hawkes_ogata`.
# #' @param t A numeric or vector.
# #'
# #' @return The compensator at time t.
# #'
# #' @export
# #'
# #' @examples
# #' # Simulate an exponential Hawkes process with baseline intensity 1 and
# #' # excitation function 1*exp(-2t)
# #' x <- hawkes_ogata(10, 1, 1, 2)
# #' plot(x)
# #' compensator_ogata(x, 0:10)
# compensator_ogata <- function(hawkes, t) {
# 	# If outside of bounds, return error
# 	if (any(t < 0) || any(t > hawkes$end))
# 		stop("t is out of bound.")
# 	# Create vector of past closest points of hawkes
# 	index <- sapply(t, function(j) {
# 		Position(function(p) p >= j, hawkes$p, nomatch=hawkes$n+1)-1
# 	})
# 	int <- hawkes$lambda * t
# 	pind <- (index > 0)
# 	int[!pind] <- int[!pind] +
# 		(hawkes$lambda0 - hawkes$lambda) / hawkes$beta * (1 - exp(-hawkes$beta * t[!pind]))
# 	int[pind] <- int[pind] + (hawkes$lambda0 - hawkes$lambda) / hawkes$beta +
# 		hawkes$alpha / hawkes$beta *
# 		(index[pind] -
# 		exp(-hawkes$beta * (t[pind] - hawkes$p[index[pind]])) * (hawkes$A[index[pind]] + 1))
# 	return( int )
# }

# #' Residuals of a simulated Hawkes process
# #'
# #' Outputs the residuals (values of the compensator at the times of arrival) of a simulated Hawkes process.
# #' Useful function for diagnosis through the random time change theorem: the residuals should follow
# #' a unit rate Poisson process.
# #'
# #' @param hawkes A simulated Hawkes process from `hawkes_ogata`.
# #'
# #' @return The intensity at time t.
# #'
# #' @export
# #'
# #' @examples
# #' # Simulate an exponential Hawkes process with baseline intensity 1 and
# #' # excitation function 1*exp(-2t)
# #' x <- hawkes_ogata(10, 1, 1, 2)
# #' plot(x)
# #' z = residuals_ogata(x)
# #' plot(z)
# #' abline(0, 1)
# residuals_ogata <- function(hawkes) {
# 	return( hawkes$lambda * hawkes$p + (hawkes$lambda0 - hawkes$lambda) / hawkes$beta +
# 		hawkes$alpha / hawkes$beta * (1:hawkes$n - 1 - hawkes$A) )
# }
