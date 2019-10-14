##############################################################################
#                       SIMULATION OF TWINSTIM PROCESS                       #
##############################################################################

# Inhomonogeous Poisson by thinning (Ogata's modified thinning algorithm)
#' @export
inhpois <- function(T, fun, M, ...) UseMethod("inhpois")
	
#' @export
inhpois.default <- function(T, fun, M, ...) {

	# Intensity function
	int <- function(t) {return( fun(t, ...) )}

	# Simulate an homonogeous Poisson on [0,T]x[0,M]
	m <- rpois(1, lambda=M*T)
	x <- runif(m, min=0, max=T)
	y <- runif(m, min=0, max=M)

	# Ogata's thinning algorithm
	accepted <- ifelse(int(x) > y, TRUE, FALSE)
	p <- sort(x[accepted])

	# Object to return
	sim <- list(p=p, T=T, fun=int, M=M, n=length(p),
			x=x, y=y, accepted=accepted, m=m, call=match.call())
	class(sim) <- "inhpois"
	return( sim )

}

#' @export
plot.inhpois <- function(inhpois, precision=1e3, ...) {
	# Conditional intensity
	matplot(z <- seq(0, inhpois$T, by=inhpois$T / precision), sapply(z, inhpois$fun), 
		type="l", ylim=c(0, max(inhpois$M)),
		xlab=expression(italic(t)), ylab=expression(italic(U)), ...)
	# Upper bound M of Ogata's modified thinning algorithm
	segments(x0=0, x1=inhpois$T, 
		y0=inhpois$M, 
		lwd=4, col="blue")
	# All points considered and the corresponding value for U
	points(x=inhpois$x, y=inhpois$y, 
		pch=ifelse(inhpois$accepted, 1, 3),
		col=ifelse(inhpois$accepted, "green4", "firebrick1"))
	# The resulting points of the Hawkes process
	points(x=inhpois$p, y=rep(0, inhpois$n), col="green4", pch=15)
	# Nice lines showing which point considered resulted in which point of the Hawkes process
	segments(x0=inhpois$x[inhpois$accepted], y0=rep(0, inhpois$n), 
		y1=inhpois$y[inhpois$accepted],
		col="green4", lty=2)
	# Legends
	legend("topleft", legend=c(expression(italic(lambda(t))), expression(italic(M)), expression(Accepted), expression(Rejected)),
		lty=c(1,1,NA,NA), col=c("black", "blue", "green", "red"), lwd=c(1,4,NA,NA), pch=c(NA, NA, 1, 3))
}

#' @export
twinstim <- function(T, fun, M, alpha, beta, ...) UseMethod("twinstim")

#' @export
twinstim.default <- function(T, fun, M, alpha, beta, ...) {
	
	# Get immigrants distributed as inhpois 
	# and number of descendants for each distributed as Pois(alpha/beta)
	immigrants <- inhpois(T, fun, M, ...)

	# Initialize
	gen <- list(gen0=immigrants$p)
	ancestors <- list()
	it <- 0

	# Generate descendants
	while (!is.null(gen[[paste0("gen", it)]])) {
		for (Ci in gen[[paste0("gen", it)]]) {
			Di <- rpois(1, lambda=alpha/beta)
			if (Di > 0) {
				Ei <- rexp(Di, rate=beta)
				gen[[paste0("gen", it+1)]] <- c(gen[[paste0("gen", it+1)]], (Ci + Ei)[(Ci + Ei) < T])
				ancestors[[paste0("gen", it+1)]] <- 
					c(ancestors[[paste0("gen", it+1)]], rep(which(gen[[paste0("gen", it)]] == Ci), Di)[(Ci + Ei) < T])
			}
		}
		it <- it + 1
	}

	# Remove last generation if empty 
	# (happens when candidates are further than T)
	if (length(gen[[paste0("gen", it-1)]]) == 0) {
		gen[[paste0("gen", it-1)]] <- NULL
		ancestors[[paste0("gen", it-1)]] <- NULL
	}
	
	# Add immigrants and sort
	if (length(gen) > 0) p <- sort(unlist(gen, use.names=FALSE))
	else p <- numeric(0)

	# Compute A
	if (length(p)== 0) A <- numeric(0)
	if (length(p) > 0) A <- 0
	if (length(p) > 1) 
		for (i in 2:length(p)) A[i] <- exp(-beta * (p[i] - p[i-1])) * (1 + A[i-1])

	# Return object
	sim <- list(p=p, T=T, alpha=alpha, beta=beta, n=length(p), A=A,
		immigrants=immigrants, gen=gen, ancestors=ancestors, call=match.call())
	class(sim) <- "twinstim"
	return( sim )

}

#' @export
plot.twinstim <- function(twinstim, intensity=FALSE, precision=1e3, ...) {
	if (intensity==FALSE) {
		# Draw a convenient empty plot
		plot(x=NULL, xlim=c(0, twinstim$T), ylim=c(-1, length(twinstim$gen)-.5), yaxt="n", 
			ylab="Generations", xlab="Time")
		axis(2, at=1:length(twinstim$gen)-1)
		# Add generation 0 (i.e. immigrants)
		points(x=twinstim$gen$gen0, y=rep(0, length(twinstim$gen$gen0)), pch=15)
		# Add further generations
		eps <- 0.05
		if (length(twinstim$gen) == 1) return()
		for (i in 2:length(twinstim$gen)) {
			points(x=twinstim$gen[[i]], y=rep(i-1, length(twinstim$gen[[i]])), pch=19, col=i)
			arrows(x0=twinstim$gen[[i]], x1=twinstim$gen[[i-1]][twinstim$ancestors[[i-1]]],
				y0=i-1-eps, y1=i-2+eps, code=1, length=.1)
		}
		# Add full process
		points(x=twinstim$p, y=rep(-0.5, twinstim$n), pch=4)
		segments(x0=0, x1=twinstim$T, y0=-.5, col="grey")
	}
	if (intensity==TRUE) {
		# Conditional intensity
		matplot(z <- seq(0, twinstim$T, by=twinstim$T / precision), 
			zt <- sapply(z, function(i) {intensity(twinstim, i)}), 
			type="l", ylim=c(0, max(zt)),
			xlab=expression(italic(t)), ylab=expression("Conditionnal intensity"), ...)
		# Hawkes process
		segments(x0=0, x1=twinstim$T, y0=0, col="grey")
		points(x=twinstim$p, y=rep(0, twinstim$n), pch=4)
	}
}

#' @export
intensity <- function(object, t) UseMethod("intensity")

#' @export
intensity.twinstim <- function(twinstim, t) {
	# If outside of bounds, return error
	if (any(t < 0) || any(t > twinstim$T)) 
		stop("t is out of bound.")
	# Create vector of past closest points of twinstim
	index <- sapply(t, function(j) {
		Position(function(p) p >= j, twinstim$p, nomatch=twinstim$n+1)-1
	})
	# Endemic part of the intensity
	int <- twinstim$immigrants$fun(t)
	# Epidemic part
	pind <- which(index > 0)
	int[pind] <- int[pind] + 
		twinstim$alpha * (twinstim$A[index[pind]]+1) * 
		exp(- twinstim$beta * (t[pind] - twinstim$p[index[pind]]))
	return( int )
}

#' @export
compensator <- function(object, t, ...) UseMethod("compensator")

#' @export
compensator.twinstim <- function(twinstim, t, ...) {
	# If outside of bounds, return error
	if (any(t < 0) || any(t > twinstim$T)) 
		stop("t is out of bound.")
	# Create vector of past closest points of twinstim
	index <- sapply(t, function(j) {
		Position(function(p) p >= j, twinstim$p, nomatch=twinstim$n+1)-1
	})
	int <- sapply(t, function(s) {
		integrate(twinstim$immigrants$fun, 0, s, ...)$value
	})
	pind <- which(index > 0)
	int[pind] <- int[pind] + 
		twinstim$alpha / twinstim$beta * 
		(index[pind] - 
		exp(-twinstim$beta * (t[pind] - twinstim$p[index[pind]])) * (twinstim$A[index[pind]] + 1))
	return( int )

}

#' @export
residuals.twinstim <- function(twinstim, ...) {
	sapply(1:twinstim$n, function(i) {
		integrate(twinstim$immigrants$fun, 0, twinstim$p[i], ...)$value + twinstim$alpha / twinstim$beta * (i - 1 - twinstim$A[i])
	})
}
