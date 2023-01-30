#' @title The power law distribution
#' @name dpowerlaw
#'
#' @description Density, distribution function, quantile function and random generation for
#' the power law distribution with shape equal to `shape` and scale equal to
#' `scale`.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param shape parameter of shape.
#' @param scale parameter of scale.
#'
#' @details
#' The density function of the power law distribution is
#' \deqn{f(t) = \theta a^\theta (a+t)^{-\theta-1}}
#' where \eqn{\theta} is the shape parameter, and \eqn{a} the scale parameter.
#'
#' @return
#' `dpowerlaw` gives the density, `ppowerlaw` gives the distribution function,
#' `qpowerlaw` gives the quantile function, and `rpowerlaw` generates random
#' deviates.
NULL

#' @rdname dpowerlaw
#' @export
dpowerlaw = function(x, shape = 1.0, scale = 1.0) {
    # shape is theta, scale is a

    if (shape <= 0.0 || scale <= 0.0)
         stop("Both shape and scale must be positive.")

    return(
        shape * exp(shape * log(scale) - (shape + 1) * log(scale + x))
    )
}

#' @rdname dpowerlaw
#' @export
ppowerlaw = function(q, shape = 1.0, scale = 1.0) {
    # shape is theta, scale is a

    if (shape <= 0.0 || scale <= 0.0)
        stop("Both shape and scale must be positive.")

    return(
        1 - exp(shape * (log(scale) - log(scale + q)))
    )
}

#' @rdname dpowerlaw
#' @export
qpowerlaw = function(p, shape = 1.0, scale = 1.0) {
    # shape is theta, scale is a

    if (shape <= 0.0 || scale <= 0.0)
        stop("Both shape and scale must be positive.")

    return(
        scale * (exp(- log(1-p)/shape) - 1)
    )
}

#' @rdname dpowerlaw
#' @export
rpowerlaw = function(n, shape = 1.0, scale = 1.0) {
    # shape is theta, scale is a

    if (shape <= 0.0 || scale <= 0.0)
        stop("Both shape and scale must be positive.")

    u = runif(n)

    return(
        scale * (exp(- log(1-u)/shape) - 1)
    )
}
