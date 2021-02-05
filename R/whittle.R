#' Fitting Hawkes processes from discrete data
#'
#' This function fits a Hawkes process to discrete data by minimizing the Whittle contrast.
#'
#' @param counts A bin-count sequence
#' @param kern Either a string (partially) matching one of the kernels implemented (see Details), or an object of class Model
#' @param binsize (Optional) The bin size of the bin-count sequence; ignored if `kern` is of class Model
#' @param init (Optional) Initial values of the optimisation algorithm
#' @param trunc (Optional) The number of foldings taken into account due to aliasing
#' @param ... Additional arguments passed to `optim`
#'
#' @details
#' If specified as string, the argument `kern` must match (partially) one of the following
#' (upper cases not taken into account): Exponential, SymmetricExponential,
#' Gaussian, PowerLaw, Pareto3, Pareto2, Pareto1.
#' The periodogram used in the optimisation procedure is computed in complexity
#' \eqn{O(n \log n)}, using function `fft`.
#'
#' @return
#' eturns a list containing the solution of the optimisation procedure, the object Model
#' with its parameters updated to the solution, and the output produced by `optim`.
#' @export
#'
#' @examples
#' # Simulate and fit a Hawkes process with exponential kernel
#' x = hawkes(100, fun = 1, repr = .5, family = "exp", rate = 1)
#' y = discrete(x, binsize = 1)
#' whittle(y, "Exponential")
#'
#' # Simulate and fit a Hawkes process with power law kernel
#' x = hawkes(100, fun = 1, repr= .3, family = "powerlaw", shape = 3.5, scale = 1.0)
#' y = discrete(x, binsize = 2)
#' whittle(y, "powerlaw", 2, init = c(2.0, .5, 3.0, 1.5), trunc = 2L)
whittle <- function(counts, kern, binsize = 1.0, init = NULL, trunc = 5L, ...) {

    # Check that the argument 'kern' is either a string that matches of the kernels implemented
    if (is.character(kern)) {
        kern = match.arg(tolower(kern),
                         c("exponential", "symmetricexponential", "gaussian",
                           "powerlaw", "pareto3", "pareto2", "pareto1"))
        switch(kern,
               exponential = {model = new(Exponential)},
               symmetricexponential = {model = new(SymmetricExponential)},
               gaussian = {model = new(Gaussian)},
               powerlaw = {model = new(PowerLaw)},
               pareto3 = {model = new(Pareto3)},
               pareto2 = {model = new(Pareto2)},
               pareto1 = {model = new(Pareto1)})
        model$binsize = binsize
    } else if ( # or that it refers to a valid hawkes model
        !any(sapply(
            paste0("Rcpp_", c("Exponential", "SymmetricExponential", "Gaussian", "PowerLaw", "Pareto3", "Pareto2", "Pareto1")),
            function(class_) {is(model, class_)}
        ))
    ) stop("'kern' must be a valid kernel.")

    # Periodogram
    n <- length(counts)
    dft <- fft(counts - mean(counts))
    I <- Mod(dft)^2 / n
    I <- I[-1]      # remove value corresponding to density in 0 (not well defined for centered processes)

    # Whittle pseudo likelihood function (for optim)
    wlik <- function(param) {
        model$param <- param
        return( model$whittle(I, trunc) )
    }

    if (is.null(init))
        x0 = c(runif(1, 0, 2),
               runif(1, 0, 1),
               runif(length(model$param)-2, 0, 2))
    else x0 = init

    optargs = list(hessian = TRUE,
                   lower = rep(.0001, length(model$param)),
                   upper = c(Inf, .9999, rep(Inf, length(model$param)-2)),
                   method = "L-BFGS-B")
    optargs = modifyList(optargs, list(...))
    opt <- do.call(optim, c(list(par = x0, fn = wlik), optargs))

    # Create output object
    model$param = opt$par

    output = list(
        par = opt$par,
        model = model,
        counts = counts,
        binsize = model$binsize,
        opt = opt
    )

    class(output) = "HawkesModel"

    return( output )
}
