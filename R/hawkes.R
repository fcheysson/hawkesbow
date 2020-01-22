##############################################################################
#                        SIMULATION OF HAWKES PROCESS                        #
##############################################################################

# Inhomonogeous Poisson by thinning (Ogata's modified thinning algorithm)

#' Simulation of an inhomogeneous Poisson process by thinning
#'
#' Simulates an inhomogeneous Poisson process via Ogata's modified thinning algorithm on [0,T].
#' An homogeneous Poisson process with intensity `M` is first generated on [0,T], then thinned using the specified intensity function `fun`.
#'
#' @param T Right bound of the interval [0,T].
#' @param fun (default = NULL) Intensity function.
#' @param M (default = NULL) Maximum bound on `fun` or, if `fun` is `NULL`, constant intensity function.
#'
#' @return A S3 object of class `inhpois` containing a vector ($p) of simulated values,
#' and all other objects used for the simulation.
#'
#' @export
#'
#' @examples
#' # Simulate an inhomogeneous Poisson process with function intensity 1 + sin(x) (bounded by 2)
#' x <- inhpois(T=10, fun=function(y) {1 + sin(y)}, M=2)
#' # Simulate a homogeneous Poisson process with intensity 3
#' x <- inhpois(T=10, M=3)
inhpois <- function(T, fun=NULL, M=NULL) UseMethod("inhpois")

#' @export
inhpois.default <- function(T, fun=NULL, M=NULL) {

    if (is.null(fun) & is.null(M)) stop("One of 'fun' or 'M' must not be 'NULL'")
    if (is.null(fun)) fun = function(x) M
    if (is.null(M)) {
        warning("M is computed approximately and should be user-specified.")
        M = max( sapply(seq(0, T, length.out = 1e3), fun) )
    }

    # Simulate an homonogeous Poisson on [0,T]x[0,M]
    m <- rpois(1, lambda=M*T)
    x <- runif(m, min=0, max=T)
    y <- runif(m, min=0, max=M)

    # Ogata's thinning algorithm
    accepted <- ifelse(fun(x) > y, TRUE, FALSE)
    p <- sort(x[accepted])

    # Object to return
    sim <- list(p=p, T=T, fun=fun, M=M, n=length(p),
                x=x, y=y, accepted=accepted, m=m, call=match.call())
    class(sim) <- "inhpois"
    return( sim )

}

#' Plot of a simulated inhomogeneous Poisson process
#'
#' Plots a simulated inhomogeneous Poisson process, highlighting the steps used in Ogata's thinning algorithm.
#'
#' @param inhpois A simulated inhomogeneous Poisson process.
#' @param precision (default = 1e3) Number of points to plot.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # Simulate an inhomogeneous Poisson process with function intensity 1 + sin(x)
#' x <- inhpois(T=10, fun=function(y) {1 + sin(y)}, M=2)
#' plot(x)
plot.inhpois <- function(inhpois, precision=1e3) {
    # Conditional intensity
    matplot(z <- seq(0, inhpois$T, length.out=precision), sapply(z, inhpois$fun),
            type="l", ylim=c(0, inhpois$M),
            xlab=expression(italic(t)), ylab=expression(italic(U)))
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
           lty=c(1,1,NA,NA), col=c("black", "blue", "green", "red"), lwd=c(1,4,NA,NA), pch=c(NA, NA, 1, 3), cex=1, horiz=FALSE)
}

#' Simulation of a Hawkes process
#'
#' Simulates a Hawkes process using its cluster representation:
#' - First, the immigrants are drawn according to an inhomogeneous Poisson process with intensity measure `fun`.
#' - Second, the number of offsprings of an immigrant follows a Poisson distribution with intensity `repr`.
#' - Third, these offsprings are distributed according to the `family` distribution.
#' - Then, generate further offsprings according to the last two steps.
#'
#' @param T A numeric value - The right bound of the interval [0,T].
#' @param fun (default = NULL) A numeric function - The intensity function of the immigrant process.
#' @param M (default = NULL) A non-negative numeric value - The maximum bound on `fun` or, if `fun` is `NULL`, constant intensity function of the immigrant process.
#' @param repr A non-negative numeric value - The mean number of offsprings.
#' @param family A character string "name" naming a distribution with corresponding random generation function `rname`, or directly the random generation function.
#' @param ... Additional arguments passed on to the random generation function `rname`.
#'
#' @return A S3 object of class Hawkes containing a vector ($p) of simulated values,
#' and all other objects used for the simulation.
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1, reproduction mean 0.5 and exponential fertility function with rate 2.
#' x <- hawkes(10, M=1, repr=0.5, family="exp", rate=2)
#' # Simulate a Hawkes process with baseline intensity function 1 + sin(x), reproduction mean 0.5 and custom [0,1]-triangular fertility function.
#' x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5, family=function(n) {1 - sqrt(1 - runif(n))})
hawkes <- function(T, fun=NULL, M=NULL, repr, family="exp", ...) UseMethod("hawkes")

#' @export
hawkes.default <- function(T, fun=NULL, M=NULL, repr, family="exp", ...) {

    # Check if fertility distribution function is user specified or chosen amongst R "r___" functions
    if (is.character(family))
        tryCatch(expr = {distfun = get(paste0("r", family), mode="function")},
                 error = function(cond) {stop(paste0("Family argument '", family, "' has no corresponding random generation function 'r", family, "'"))})
    else if (!is.function(family)) stop(paste0("The argument 'family' is not a function."))
    else distfun = family

    # Get immigrants distributed as inhpois
    # and number of descendants for each distributed as Poisson with rate `repr`
    immigrants <- inhpois(T, fun, M)

    # Initialize
    gen <- list(gen0=immigrants$p)
    ancestors <- list()
    it <- 0

    # Generate descendants
    while (!is.null(gen[[paste0("gen", it)]])) {
        for (Ci in gen[[paste0("gen", it)]]) {
            Di <- rpois(1, lambda=repr)
            if (Di > 0) {
                Ei <- distfun(Di, ...)
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

    # Return object
    sim <- list(p=p, T=T, repr=repr, distfun=distfun, family=family, distargs=list(...), n=length(p),
                immigrants=immigrants, gen=gen, ancestors=ancestors, call=match.call())
    class(sim) <- "hawkes"
    return( sim )

}

#### TO FINISHHHHHHH
#' @export
plot.hawkes <- function(hawkes, intensity=FALSE, precision=1e3) {
    if (intensity==FALSE) {
        # Draw a convenient empty plot
        plot(x=NULL, xlim=c(0, hawkes$T), ylim=c(-1, length(hawkes$gen)-.5), yaxt="n",
             ylab="Generations", xlab="Time")
        axis(2, at=1:length(hawkes$gen)-1)
        # Add generation 0 (i.e. immigrants)
        points(x=hawkes$gen$gen0, y=rep(0, length(hawkes$gen$gen0)), pch=15)
        # Add further generations
        eps <- 0.05
        if (length(hawkes$gen) == 1) return()
        for (i in 2:length(hawkes$gen)) {
            points(x=hawkes$gen[[i]], y=rep(i-1, length(hawkes$gen[[i]])), pch=19, col=i)
            arrows(x0=hawkes$gen[[i]], x1=hawkes$gen[[i-1]][hawkes$ancestors[[i-1]]],
                   y0=i-1-eps, y1=i-2+eps, code=1, length=.1)
        }
        # Add full process
        points(x=hawkes$p, y=rep(-0.5, hawkes$n), pch=4)
        segments(x0=0, x1=hawkes$T, y0=-.5, col="grey")
    }
    if (intensity==TRUE) {
        # Conditional intensity
        matplot(z <- seq(0, hawkes$T, by=hawkes$T / precision),
                zt <- sapply(z, function(i) {intensity(hawkes, i)}),
                type="l", ylim=c(0, max(zt)),
                xlab=expression(italic(t)), ylab=expression("Conditionnal intensity"))
        # Hawkes process
        segments(x0=0, x1=hawkes$T, y0=0, col="grey")
        points(x=hawkes$p, y=rep(0, hawkes$n), pch=4)
    }
}

#' Intensity of a Hawkes process
#'
#' Outputs the intensity of a Hawkes process `x`, given the specified set of parameters.
#' If the input `x` is of class "hawkes" output by function `hawkes`, the parameters of the simulation will be used to compute the intensity.
#' If any parameter is user-specified in the function call, it will use these instead.
#'
#' @param x Either: a numeric vector, sorted in ascending order; or an object of class "hawkes" output by function `hawkes`.
#' @param t A numeric value or vector, at which the intensity should be computed.
#' @param fun (default = NULL) A numeric function - The intensity function of the immigrant process.
#' @param M (default = NULL) A non-negative numeric value - The maximum bound on `fun` or, if `fun` is `NULL`, constant intensity function of the immigrant process.
#' @param repr (default = NULL) A non-negative numeric value - The mean number of offsprings.
#' @param family (default = NULL) A character string "name" naming a distribution with corresponding distribution function `dname`, or directly the distribution function.
#' @param ... Additional arguments passed on to the random generation function `dname`.
#'
#' @return The intensity at time t.
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1, reproduction mean 0.5 and exponential fertility distribution with rate 2.
#' x <- hawkes(10, M=1, repr=0.5, family="exp", rate=2)
#' plot(x, intensity = TRUE)
#' intensity(x, 0:10)
#' # Intensity under a different model
#' intensity(x, 0:10, repr=0.8, rate=3)
#' # Simulate a Hawkes process with baseline intensity function 1 + sin(x), reproduction mean 0.5 and custom [0,1]-triangular fertility function.
#' x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5, family=function(n) {1 - sqrt(1 - runif(n))})
#' plot(x, intensity = TRUE)
#' intensity(x, 0:10, family=function(y) ifelse(y>0 & y<1, 2-2*y, 0))
intensity <- function(x, t, fun = NULL, M = NULL, repr = NULL, family = NULL, ...) {

    # Define pattern
    if (is.numeric(x)) p = x
    else if (is(x, "hawkes")) p = x$p
    else stop("'x' must be numeric.")

    # Define immigrant intensity function
    if (!is.null(fun)) {}
    else if (!is.null(M)) fun = function(y) M
    else if (is(x, "hawkes")) fun = x$immigrants$fun
    else stop("One of 'fun' or 'M' must not be 'NULL'.")

    # Define reproduction mean
    if (!is.null(repr)) {}
    else if (is(x, "hawkes")) repr = x$repr
    else stop("'repr' must not be 'NULL'.")

    # Define family
    if (is.character(family))
        tryCatch(expr = {distfun = get(paste0("d", family), mode="function")},
                 error = function(cond) {stop("d", family, " is not a valid density function.")})
    else if (!is.null(family) & !is.function(family)) stop(paste0("The argument 'family' is not a function."))
    else if (!is.null(family)) {
        distfun = family
        family = deparse(substitute(family))
    } else if (is(x, "hawkes")) {
        family = x$family
        tryCatch(expr = {distfun = get(paste0("d", x$family), mode="function")},
                 error = function(cond) {stop("d", x$family, " is not a valid density function. Specify one in argument 'family'.")})
    } else stop("'family' must not be 'NULL'.")

    # Define optional arguments if user-specified (or overwrite if 'x' is of class 'hawkes')
    if (is(x, "hawkes")) distargs = modifyList(x$distargs, list(...))
    else distargs = list(...)

    # If outside of bounds, return error
    if (any(t < 0))
        stop("t must be non negative.")

    if (family == "exp") {

        rate = ifelse("rate" %in% names(distargs), distargs$rate, 1)

        # Compute A
        if (length(p)== 0) A <- numeric(0)
        if (length(p) > 0) A <- 0
        if (length(p) > 1)
            for (i in 2:length(p)) A[i] <- exp(-rate * (p[i] - p[i-1])) * (1 + A[i-1])

        # Create vector of past closest points of hawkes
        index <- sapply(t, function(j) {
            Position(function(u) u >= j, p, nomatch=length(p)+1)-1
        })

        # Endemic part of the intensity
        int <- sapply(t, fun)

        # Epidemic part
        pind <- which(index > 0)
        int[pind] <- int[pind] + repr * rate * (A[index[pind]]+1) * exp(- rate * (t[pind] - p[index[pind]]))

        return(int)

    } else {

        # Create vector of past closest points of hawkes
        index <- sapply(t, function(j) {
            Position(function(u) u >= j, p, nomatch=length(p)+1)-1
        })

        # Endemic part of the intensity
        int <- sapply(t, fun)

        # Epidemic part
        pind <- which(index > 0)
        int[pind] <- int[pind] + repr * sapply(pind, function(s) {
            sum(sapply( p[1:index[s]], function(q) {do.call(distfun, c(t[s]-q, distargs))} ))
        })

        return(int)

    }

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
