##############################################################################
#                        SIMULATION OF HAWKES PROCESS                        #
##############################################################################

# Inhomonogeous Poisson by thinning (Ogata's modified thinning algorithm)

#' Simulation of an inhomogeneous Poisson process by thinning
#'
#' Simulates an inhomogeneous Poisson process via Ogata's modified thinning algorithm on \eqn{[0,\mathrm{end}]}.
#' An homogeneous Poisson process with intensity `M` is first generated on \eqn{[0,\mathrm{end}]}, then thinned using the specified intensity function `fun`.
#'
#' @param end A non-negative numeric value - right bound of the interval \eqn{[0,\mathrm{end}]}.
#' @param fun A non-negative function or numeric value - intensity (function) of the Poisson process.
#' @param M (default = NULL) A non-negative numeric value - upper bound on `fun` (ignored if `fun` is a numeric value).
#'
#' @return A S3 object of class `inhpois` containing a vector ($p) of simulated values,
#' and all other objects used for the simulation.
#'
#' @export
#'
#' @examples
#' # Simulate an inhomogeneous Poisson process with function intensity 1 + sin(x) (bounded by 2)
#' x <- inhpois(end=10, fun=function(y) {1 + sin(y)}, M=2)
#' # Simulate a homogeneous Poisson process with intensity 3
#' x <- inhpois(end=10, fun=3)
inhpois <- function(end, fun, M=NULL) {

    if (is.numeric(fun)) {
        m = rpois(1, lambda = fun*end)
        p <- sort(runif(m, min = 0, max = end))

        # Object to return
        sim = list(p = p, end = end, fun = fun, n = m, call=match.call())
        return(sim)
    }

    if (!is.function(fun)) stop("Parameter 'fun' must be a non-negative function or numeric.")

    if (is.null(M)) {
        warning("M is computed approximately and should be user-specified.")
        M = max( sapply(seq(0, end, length.out = 1e3), fun) )
    }

    # Simulate an homonogeous Poisson on [0,end]x[0,M]
    m <- rpois(1, lambda=M*end)
    x <- runif(m, min=0, max=end)
    y <- runif(m, min=0, max=M)

    # Ogata's thinning algorithm
    accepted <- ifelse(fun(x) > y, TRUE, FALSE)
    p <- sort(x[accepted])

    # Object to return
    sim <- list(p=p, end=end, fun=fun, M=M, n=length(p),
                x=x, y=y, accepted=accepted, m=m, call=match.call())
    class(sim) <- "inhpois"
    return(sim)

}

#' Plot of a simulated inhomogeneous Poisson process
#'
#' Plots a simulated inhomogeneous Poisson process, highlighting the steps used in Ogata's thinning algorithm.
#'
#' @param x A simulated inhomogeneous Poisson process.
#' @param precision (default = 1e3) Number of points to plot.
#' @param ... Only there to fit the declaration of S3 method `plot`.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # Simulate an inhomogeneous Poisson process with function intensity 1 + sin(x)
#' x <- inhpois(end=10, fun=function(y) {1 + sin(y)}, M=2)
#' plot(x)
plot.inhpois <- function(x, precision=1e3, ...) {
    # Conditional intensity
    matplot(z <- seq(0, x$end, length.out=precision), sapply(z, x$fun),
            type="l", ylim=c(0, x$M),
            xlab=expression(italic(t)), ylab=expression(italic(U)))
    # Upper bound M of Ogata's modified thinning algorithm
    segments(x0=0, x1=x$end,
             y0=x$M,
             lwd=4, col="blue")
    # All points considered and the corresponding value for U
    points(x=x$x, y=x$y,
           pch=ifelse(x$accepted, 1, 3),
           col=ifelse(x$accepted, "green4", "firebrick1"))
    # The resulting points of the Hawkes process
    points(x=x$p, y=rep(0, x$n), col="green4", pch=15)
    # Nice lines showing which point considered resulted in which point of the Hawkes process
    segments(x0=x$x[x$accepted], y0=rep(0, x$n),
             y1=x$y[x$accepted],
             col="green4", lty=2)
    # Legends
    legend("topleft", legend=c(expression(italic(lambda(t))), expression(italic(M)), expression(Accepted), expression(Rejected)),
           lty=c(1,1,NA,NA), col=c("black", "blue", "green", "red"), lwd=c(1,4,NA,NA), pch=c(NA, NA, 1, 3), cex=1, horiz=FALSE)
}

#' Simulation of a Hawkes process
#'
#' Simulates a Hawkes process using its cluster representation:
#' - First, the immigrants are drawn according to an (inhomogeneous) Poisson process with intensity measure `fun`.
#' - Second, the number of offsprings of an immigrant follows a Poisson distribution with intensity `repr`.
#' - Third, these offsprings are distributed according to the `family` distribution.
#' - Then, generate further offsprings according to the last two steps.
#'
#' @param end A non-negative numeric value - right bound of the interval \eqn{[0,\mathrm{end}]}.
#' @param fun A non-negative function or numeric value - intensity (function) of the immigrant process.
#' @param repr A non-negative numeric value - mean number of offsprings.
#' @param family A character string "name" naming a distribution with corresponding random generation function `rname`, or directly the random generation function.
#' @param M (default = NULL) A non-negative numeric value - upper bound on `fun`(ignored if `fun` is a numeric value).
#' @param ... Additional arguments passed on to the random generation function.
#'
#' @return A S3 object of class Hawkes containing a vector ($p) of simulated values,
#' and all other objects used for the simulation.
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1,
#' # reproduction mean 0.5 and exponential fertility function with rate 2.
#' x <- hawkes(10, fun=1, repr=0.5, family="exp", rate=2)
#' # Simulate a Hawkes process with baseline intensity function 1 + sin(x),
#' # reproduction mean 0.5 and custom [0,1]-triangular fertility function.
#' x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5,
#'             family=function(n) {1 - sqrt(1 - runif(n))})
hawkes <- function(end, fun, repr, family, M=NULL, ...) {

    # Check if fertility distribution function is user specified or chosen amongst R "r___" functions
    if (is.character(family))
        tryCatch(expr = {distfun = get(paste0("r", family), mode="function")},
                 error = function(cond) {stop(paste0("Family argument '", family, "' has no corresponding random generation function 'r", family, "'"))})
    else if (!is.function(family)) stop(paste0("The argument 'family' is not a function."))
    else distfun = family

    # Get immigrants distributed as inhpois
    immigrants <- inhpois(end, fun, M)

    # Initialize
    gen <- list(gen0=immigrants$p)
    ancestors <- list()
    it <- 0

    # Generate offsprings
    while (length(gen[[paste0("gen", it)]])) {
        Ci <- gen[[paste0("gen", it)]]
        Di <- rpois(length(Ci), lambda=repr)
        Ei <- lapply(Di, distfun, ...)
        Fi <- lapply(1:length(Ci), function(i) {
            offsprings <- Ci[i] + Ei[[i]]
            offsprings
        })
        gen[[paste0("gen", it+1)]] <- unlist(Fi)
        ancestors[[paste0("gen", it+1)]] <- rep.int(1:length(Ci), times=sapply(Fi, length))
        it <- it + 1
    }

    # Remove last generation that is always empty (stopping criterion for while loop)
    if (length(gen[[paste0("gen", it)]]) == 0) {
        gen[[paste0("gen", it)]] <- NULL
        ancestors[[paste0("gen", it)]] <- NULL
    }

    # Add immigrants and sort
    if (length(gen) > 0) {
        p <- sort(unlist(gen, use.names=FALSE))
        p <- p[p < end & p > 0]
    } else p <- numeric(0)

    # Return object
    sim <- list(p=p, end=end, repr=repr, distfun=distfun, family=family, distargs=list(...), n=length(p),
                immigrants=immigrants, gen=gen, ancestors=ancestors, call=match.call())
    class(sim) <- "hawkes"
    return( sim )

}

#' Plot of a Hawkes process
#'
#' Plots the realisation of a Hawkes process and either its cluster representation (`intensity=FALSE`, only available for a simulated Hawkes process) or its intensity function (`intensity=TRUE`).
#'
#' @param x Either: a numeric vector, sorted in ascending order; or an object of class "hawkes" output by function `hawkes`.
#' @param intensity (default = FALSE) A boolean - whether to represent the cluster representation (`FALSE`) or the intensity function (`TRUE`).
#' @param precision (default = 1e3) Number of points to plot.
#' @param fun (default = NULL) A numeric function - intensity (function) of the immigrant process.
#' @param repr (default = NULL) A non-negative numeric value - mean number of offsprings.
#' @param family (default = NULL) A character string "name" naming a distribution with corresponding distribution function `dname`, or directly the distribution function.
#' @param M (default = NULL) A non-negative numeric value - upper bound on `fun`(ignored if `fun` is a numeric value).
#' @param ... Additional arguments passed on to the random generation function `dname`.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1,
#' # reproduction mean 0.5 and exponential fertility function with rate 2.
#' x <- hawkes(10, fun=1, repr=0.5, family="exp", rate=2)
#' plot(x)
#' # Simulate a Hawkes process with baseline intensity function 1 + sin(x),
#' # reproduction mean 0.5 and custom [0,1]-triangular fertility function.
#' x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5,
#'             family=function(n) {1 - sqrt(1 - runif(n))})
#' plot(x, intensity=TRUE, family=function(y) ifelse(y>0 & y<1, 2-2*y, 0))
plot.hawkes <- function(x, intensity = FALSE, precision = 1e3, fun = NULL, repr = NULL, family = NULL, M = NULL, ...) {
    if (intensity==FALSE) {
        if (!is(x, "hawkes")) stop("'intensity==FALSE' is only compatible with Hawkes processes simulated from function 'hawkes'.")
        # Draw a convenient empty plot
        plot(x=NULL, xlim=c(0, x$end), ylim=c(-1, length(x$gen)-.5), yaxt="n",
             ylab="Generations", xlab="Time")
        axis(2, at=1:length(x$gen)-1)
        # Add generation 0 (i.e. immigrants)
        points(x=x$gen$gen0, y=rep(0, length(x$gen$gen0)), pch=15)
        # Add further generations
        eps <- 0.05
        if (length(x$gen) == 1) return()
        for (i in 2:length(x$gen)) {
            points(x=x$gen[[i]], y=rep(i-1, length(x$gen[[i]])), pch=19, col=i)
            arrows(x0=x$gen[[i]], x1=x$gen[[i-1]][x$ancestors[[i-1]]],
                   y0=i-1-eps, y1=i-2+eps, code=1, length=.1)
        }
        # Add full process
        points(x=x$p, y=rep(-0.5, x$n), pch=4)
        segments(x0=0, x1=x$end, y0=-.5, col="grey")
    }
    if (intensity==TRUE) {
        # Define pattern
        if (is.numeric(x)) p = x
        else if (is(x, "hawkes")) p = x$p
        else stop("Argument 'x' must either be of class 'hawkes' or a numeric vector.")

        # Define end point
        if (is(x, "hawkes")) end = x$end
        else end = tail(x, 1)

        # Conditional intensity
        matplot(z <- seq(0, end, by=end / precision),
                zt <- intensity(x = x, t = z, fun = fun, M = M, repr = repr, family = family, ...),
                type="l", ylim=c(0, max(zt)),
                xlab=expression(italic(t)), ylab=expression("Conditional intensity"))
        # Hawkes process
        segments(x0=0, x1=end, y0=0, col="grey")
        points(x=p, y=rep(0, length(p)), pch=4)
    }
}

#' Intensity of a Hawkes process
#'
#' Outputs the intensity of a Hawkes process `x`, given the specified set of parameters.
#'
#' If the input `x` has been simulated using the function `hawkes`, the parameters of the simulation will be used by default to compute the intensity.
#' If any parameter is specified in this function call, the function will use this instead.
#'
#' @param x A non-negative numeric vector, sorted in ascending order; or an object of class "hawkes" output by function `hawkes`.
#' @param t A non-negative numeric value or vector, at which the intensity should be computed.
#' @param fun (default = TRUE) A non-negative numeric function or value - intensity (function) of the immigrant process.
#' @param repr (default = NULL) A non-negative numeric value - mean number of offsprings.
#' @param family (default = NULL) A character string "name" naming a distribution with corresponding distribution function `dname`, or directly the distribution function.
#' @param M (default = NULL) A non-negative numeric value - upper bound on `fun`(ignored if `fun` is a numeric value).
#' @param ... Additional arguments passed on to the random generation function `dname`.
#'
#' @return The intensity at time t.
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1,
#' # reproduction mean 0.5 and exponential fertility distribution with rate 2.
#' x <- hawkes(10, fun=1, repr=0.5, family="exp", rate=2)
#' intensity(x, 0:10)
#' # Intensity with a different set of parameters
#' intensity(x, 0:10, repr=0.8, rate=3)
#' # Intensity with a different distribution function
#' intensity(x, 0:10, family="chisq", df=2)
#' # Simulate a Hawkes process with baseline intensity function 1 + sin(x),
#' # reproduction mean 0.5 and custom [0,1]-triangular fertility function.
#' x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5,
#'             family=function(n) {1 - sqrt(1 - runif(n))})
#' intensity(x, 0:10, family=function(y) ifelse(y>0 & y<1, 2-2*y, 0))
intensity <- function(x, t, fun = NULL, repr = NULL, family = NULL, M = NULL, ...) {

    # Define pattern
    if (is.numeric(x)) p = x
    else if (is(x, "hawkes")) p = x$p
    else stop("'x' must be numeric.")

    # Define immigrant intensity function
    if (!is.null(fun)) {}
    else if (is(x, "hawkes")) fun = x$immigrants$fun
    else stop("'fun' must not be 'NULL'.")

    # Define reproduction mean
    if (!is.null(repr)) {}
    else if (is(x, "hawkes")) repr = x$repr
    else stop("'repr' must not be 'NULL'.")

    # Define family
    if (is.character(family))
        tryCatch(expr = {distfun = get(paste0("d", family), mode="function")},
                 error = function(cond) {stop("d", family, " is not a valid density function.")})
    else if (!is.null(family) && !is.function(family)) stop(paste0("The argument 'family' is not a function."))
    else if (!is.null(family)) {
        distfun = family
    } else if (is(x, "hawkes")) {
        if (is.character(x$family)) {
            family = x$family
            tryCatch(expr = {distfun = get(paste0("d", x$family), mode="function")},
                     error = function(cond) {stop("The random generation function of 'x', 'r", x$family, "', cannot be matched with its associated density function, 'd", x$family, "'. Specify one in argument 'family'.")})
        } else stop("The random generation function of 'x' cannot be matched with its associated density function. Specify one in argument 'family'.")
    } else stop("'family' must not be 'NULL'.")

    # Define optional arguments if user-specified (or overwrite if 'x' is of class 'hawkes')
    if (is(x, "hawkes")) distargs = modifyList(x$distargs, list(...))
    else distargs = list(...)

    # Remove unused arguments to the distfun function
    distargs = distargs[names(distargs) %in% names(formals(distfun))]

    # If outside of bounds, return error
    if (any(t < 0))
        stop("t must be non negative.")

    if (is.character(family) && family == "exp") {

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
        if (is.numeric(fun))
            int <- rep(fun, length(t))
        else if (is.function(fun))
            int <- sapply(t, fun)
        else stop("'fun' must be a function or numeric value.")

        # Epidemic part
        pind <- which(index > 0)
        if (length(pind) > 0)
            int[pind] <- int[pind] + repr * rate * (A[index[pind]]+1) * exp(- rate * (t[pind] - p[index[pind]]))

        return(int)

    } else {

        # Create vector of past closest points of hawkes
        index <- sapply(t, function(j) {
            Position(function(u) u >= j, p, nomatch=length(p)+1)-1
        })

        # Endemic part of the intensity
        if (is.numeric(fun))
            int <- rep(fun, length(t))
        else if (is.function(fun))
            int <- sapply(t, fun)
        else stop("'fun' must be a function or numeric value.")

        # Epidemic part
        pind <- which(index > 0)
        if (length(pind) > 0) {
            int[pind] <- int[pind] + repr * sapply(pind, function(s) {
                sum(sapply( p[1:index[s]], function(q) {do.call(distfun, c(t[s]-q, distargs))} ))
            })
        }

        return(int)

    }

}

#' Compensator of a Hawkes process
#'
#' Outputs the compensator (integrated intensity) of a Hawkes process.
#'
#' @param x A non-negative numeric vector, sorted in ascending order; or an object of class "hawkes" output by function `hawkes`.
#' @param t A non-negative numeric value or vector, at which the intensity should be computed.
#' @param fun (default = TRUE) A non-negative numeric function or value - intensity (function) of the immigrant process.
#' @param repr (default = NULL) A non-negative numeric value - mean number of offsprings.
#' @param family (default = NULL) A character string "name" naming a distribution with corresponding distribution function `dname`, or directly the distribution function.
#' @param M (default = NULL) A non-negative numeric value - upper bound on `fun`(ignored if `fun` is a numeric value).
#' @param ... Additional arguments passed on to the random generation function `dname`.
#'
#' @return The compensator at time t.
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1,
#' # reproduction mean 0.5 and exponential fertility distribution with rate 2.
#' x <- hawkes(10, fun=1, repr=0.5, family="exp", rate=2)
#' compensator(x, 0:10)
#' # Compensator with a different set of parameters
#' compensator(x, 0:10, repr=0.8, rate=3)
#' # Compensator with a different distribution function
#' compensator(x, 0:10, family="chisq", df=2)
#' # Simulate a Hawkes process with baseline intensity function 1 + sin(x),
#' # reproduction mean 0.5 and custom [0,1]-triangular fertility function.
#' x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5,
#'             family=function(n) {1 - sqrt(1 - runif(n))})
#' compensator(x, 0:10, family=function(y) ifelse(y>0 & y<1, 2-2*y, 0))
compensator <- function(x, t, fun = NULL, repr = NULL, family = NULL, M = NULL, ...) {

    # Check if family should be inherited from x
    if (!is.null(family)) family_ = family
    else if (is(x, "hawkes")) family_ = x$family

    # If outside of bounds, return error
    if (any(t < 0))
        stop("t must be non negative.")

    if (is.character(family_) && family_ == "exp") {

        # Define pattern
        if (is.numeric(x)) p = x
        else if (is(x, "hawkes")) p = x$p
        else stop("'x' must be numeric.")

        # Define immigrant intensity function
        if (!is.null(fun)) {}
        else if (is(x, "hawkes")) fun = x$immigrants$fun
        else stop("'fun' must not be 'NULL'.")

        # Define reproduction mean
        if (!is.null(repr)) {}
        else if (is(x, "hawkes")) repr = x$repr
        else stop("'repr' must not be 'NULL'.")

        # Define optional arguments if user-specified (or overwrite if 'x' is of class 'hawkes')
        if (is(x, "hawkes")) distargs = modifyList(x$distargs, list(...))
        else distargs = list(...)

        rate = ifelse("rate" %in% names(distargs), distargs$rate, 1)

        # Compute A
        if (length(p)== 0) A <- numeric(0)
        if (length(p) > 0) A <- 0
        if (length(p) > 1)
            for (i in 2:length(p)) A[i] <- exp(-rate * (p[i] - p[i-1])) * (1 + A[i-1])

        # Create vector of past closest points of twinstim
        index <- sapply(t, function(j) {
            Position(function(u) u >= j, p, nomatch=length(p)+1)-1
        })

        # Endemic part of the compensator
        if (is.numeric(fun))
            int <- fun * t
        else if (is.function(fun)) {
            partition = c(0, t)
            int <- cumsum(sapply(2:length(partition), function(s) {
                integrate(function(y) sapply(y, fun), partition[s-1], partition[s])$value
            }))
        } else stop("'fun' must be a function or numeric value.")

        # Epidemic part
        pind <- which(index > 0)
        if (length(pind) > 0)
            int[pind] <- int[pind] + repr * (index[pind] - exp(- rate * (t[pind] - p[index[pind]])) * (A[index[pind]]+1))

        return( int )

    } else {

        partition = c(0, t)
        int <- cumsum(sapply(2:length(partition), function(s) {
            integrate(intensity, partition[s-1], partition[s], x=x, fun=fun, M=M, repr=repr, family=family, ...)$value
        }))

        return(int)

    }

}

#' Residuals of a Hawkes process
#'
#' Outputs the residuals (values of the compensator at the times of arrival) of a Hawkes process.
#' Useful function for diagnosis through the random time change theorem: the residuals should follow
#' a unit rate Poisson process.
#'
#' @param x A non-negative numeric vector, sorted in ascending order; or an object of class "hawkes" output by function `hawkes`.
#' @param fun (default = TRUE) A non-negative numeric function or value - intensity (function) of the immigrant process.
#' @param repr (default = NULL) A non-negative numeric value - mean number of offsprings.
#' @param family (default = NULL) A character string "name" naming a distribution with corresponding distribution function `dname`, or directly the distribution function.
#' @param M (default = NULL) A non-negative numeric value - upper bound on `fun`(ignored if `fun` is a numeric value).
#' @param ... Additional arguments passed on to the random generation function `dname`.
#'
#' @return The residuals of the process.
#'
#' @export
#'
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1,
#' # reproduction mean 0.5 and exponential fertility distribution with rate 2.
#' x <- hawkes(10, fun=1, repr=0.5, family="exp", rate=2)
#' resid = residuals(x)
#' resid
#' plot(resid)
#' abline(0, 1, col="red", lty="dashed")
#' # Residuals with a different set of parameters
#' residuals(x, repr=0.8, rate=3)
#' # Residuals with a different distribution function
#' residuals(x, family="chisq", df=2)
#' # Simulate a Hawkes process with baseline intensity function 1 + sin(x),
#' # reproduction mean 0.5 and custom [0,1]-triangular fertility function.
#' x <- hawkes(10, fun=function(y) {1+sin(y)}, M=2, repr=0.5,
#'             family=function(n) {1 - sqrt(1 - runif(n))})
#' resid = residuals(x, family=function(y) ifelse(y>0 & y<1, 2-2*y, 0))
#' plot(resid)
#' abline(0, 1, col="red", lty="dashed")
residuals <- function(x, fun = NULL, repr = NULL, family = NULL, M = NULL, ...) {

    # Check if family should be inherited from x
    if (!is.null(family)) family_ = family
    else if (is(x, "hawkes")) family_ = x$family

    # Define pattern
    if (is.numeric(x)) p = x
    else if (is(x, "hawkes")) p = x$p
    else stop("'x' must be numeric.")

    if (is.character(family_) && family_ == "exp") {

        # Define immigrant intensity function
        if (!is.null(fun)) {}
        else if (is(x, "hawkes")) fun = x$immigrants$fun
        else stop("'fun' must not be 'NULL'.")

        # Define reproduction mean
        if (!is.null(repr)) {}
        else if (is(x, "hawkes")) repr = x$repr
        else stop("'repr' must not be 'NULL'.")

        # Define optional arguments if user-specified (or overwrite if 'x' is of class 'hawkes')
        if (is(x, "hawkes")) distargs = modifyList(x$distargs, list(...))
        else distargs = list(...)

        rate = ifelse("rate" %in% names(distargs), distargs$rate, 1)

        # Compute A
        if (length(p)== 0) A <- numeric(0)
        if (length(p) > 0) A <- 0
        if (length(p) > 1)
            for (i in 2:length(p)) A[i] <- exp(-rate * (p[i] - p[i-1])) * (1 + A[i-1])

        # Endemic part of the compensator
        if (is.numeric(fun))
            resid <- fun * p
        else if (is.function(fun)) {
            partition = c(0, p)
            resid <- cumsum(sapply(2:length(partition), function(s) {
                integrate(function(y) sapply(y, fun), partition[s-1], partition[s])$value
            }))
        } else stop("'fun' must be a function or numeric value.")

        # Epidemic part
        resid = resid + sapply(1:length(p), function(i) {
            repr * (i - 1 - A[i])
        })

        return(resid)

    } else {

        resid = compensator(x = x, t = p, fun = fun, M = M, repr = repr, family = family, ...)
        return(resid)

    }

}
