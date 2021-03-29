#' Fitting Hawkes processes from discrete data
#'
#' This function fits a Hawkes process to discrete data by minimizing the Whittle contrast.
#'
#' @param counts A bin-count sequence
#' @param kern Either a string (partially) matching one of the kernels implemented (see Details), or an object of class Model
#' @param binsize (Optional) The bin size of the bin-count sequence; if omitted, defaults to 1 if `kern` is a string, or uses the member `binsize` of `kern` if it is of class Model
#' @param trunc (Optional) The number of foldings taken into account due to aliasing
#' @param init (Optional) Initial values of the optimisation algorithm
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
#' Returns a list containing the solution of the optimisation procedure, the object Model
#' with its parameters updated to the solution, and the output produced by `optim`.
#' @export
#'
#' @seealso [hawkes()] for the simulation of Hawkes processes,
#' [discrete()] for the discretisation of simulated Hawkes processes,
#' [Model] for the abstract class, and [Exponential] for the specific
#' reproduction kernels.
#'
#' @examples
#' # Simulate and fit a Hawkes process with exponential kernel
#' x = hawkes(1000, fun = 1, repr = .5, family = "exp", rate = 1)
#' y = discrete(x, binsize = 1)
#' opt = whittle(y, "Exponential")
#' opt$par      # Estimated parameters
#'
#' \donttest{
#' # May take up to 20 seconds
#' # Simulate and fit a Hawkes process with power law kernel
#' x = hawkes(1000, fun = 1, repr= .3, family = "powerlaw", shape = 3.5, scale = 1.0)
#' y = discrete(x, binsize = 1)
#' opt = whittle(y, "powerlaw")
#' opt$par      # Estimated parameters
#' }
whittle <- function(counts, kern, binsize = NULL, trunc = 5L, init = NULL, ...) {

    # Check that the argument 'kern' is either a string that matches of the kernels implemented
    if (is.character(kern)) {
        kern = match.arg(tolower(kern),
                         c("exponential", "symmetricexponential", "gaussian",
                           "powerlaw", "pareto3", "pareto2", "pareto1"))
        switch(kern,
               exponential = {kern = new(Exponential)},
               symmetricexponential = {kern = new(SymmetricExponential)},
               gaussian = {kern = new(Gaussian)},
               powerlaw = {kern = new(PowerLaw)},
               pareto3 = {kern = new(Pareto3)},
               pareto2 = {kern = new(Pareto2)},
               pareto1 = {kern = new(Pareto1)})
    } else if ( # or that it refers to a valid hawkes kernel
        !any(sapply(
            paste0("Rcpp_", c("Exponential", "SymmetricExponential", "Gaussian", "PowerLaw", "Pareto3", "Pareto2", "Pareto1")),
            function(class_) {is(kern, class_)}
        ))
    ) stop("'kern' must be a valid kernel.")

    if (!is.null(binsize)) kern$binsize = binsize

    # Periodogram
    n <- length(counts)
    dft <- fft(counts - mean(counts))
    I <- Mod(dft)^2 / n
    I <- I[-1]      # remove value corresponding to density in 0 (not well defined for centered processes)

    # Whittle pseudo likelihood function (for optim)
    nlopt_fn <- function(param) {
        kern$param <- param
        return( kern$whittle(I, trunc) )
    }

    # Sensible initialisation
    if (is.null(init)) {
        ymean = mean(counts)
        # For PowerLaw
        if (is(kern, "Rcpp_PowerLaw")) {
            wmin = Inf
            mu = .25
            shape = 1
            scale = 1
            for (mu_ in c(.25, .5, .75)) {
                for (shape_ in 1:10/2) {
                    for (scale_ in 1:10/2) {
                        kern$param[1] = ymean * (1 - mu_)
                        kern$param[2] = mu_
                        kern$param[3] = shape_
                        kern$param[4] = scale_
                        whit = kern$whittle(I, trunc = trunc)
                        if (whit < wmin) {
                            mu = mu_
                            shape = shape_
                            scale = scale_
                            wmin = whit
                        }
                    }
                }
            }
            x0 = c(ymean * (1 - mu), mu, shape, scale)
        } else if (is(kern, "Rcpp_Exponential")) {
            wmin = Inf
            mu = .25
            rate = 1
            for (mu_ in c(.25, .5, .75)) {
                for (rate_ in 1:5) {
                    kern$param[1] = ymean * (1 - mu_)
                    kern$param[2] = mu_
                    kern$param[3] = rate_
                    whit = kern$whittle(I, trunc = trunc)
                    if (whit < wmin) {
                        mu = mu_
                        rate = rate_
                        wmin = whit
                    }
                }
            }
            x0 = c(ymean * (1 - mu), mu, rate)
        } else {
            x0 = c(ymean * 2,
                   .5,
                   runif(length(kern$param)-2, 0, 10))
        }
    } else x0 = init

    optargs = list(hessian = TRUE,
                   lower = rep(.01, length(kern$param)),
                   upper = c(Inf, .99, rep(Inf, length(kern$param)-2)),
                   method = "L-BFGS-B")

    optargs = modifyList(optargs, list(...))
    opt <- do.call(optim, c(list(par = x0, fn = nlopt_fn), optargs))

    # Create output object
    kern$param = opt$par

    output = list(
        par = opt$par,
        kernel = kern,
        counts = counts,
        binsize = kern$binsize,
        opt = opt
    )

    class(output) = "HawkesModel"

    return( output )
}

# whittle2 <- function(counts, kern, binsize = NULL, trunc = 5L, init = NULL, opts = NULL, ...) {
#
#     # Check that the argument 'kern' is either a string that matches of the kernels implemented
#     if (is.character(kern)) {
#         kern = match.arg(tolower(kern),
#                          c("exponential", "symmetricexponential", "gaussian",
#                            "powerlaw", "pareto3", "pareto2", "pareto1"))
#         switch(kern,
#                exponential = {kern = new(Exponential)},
#                symmetricexponential = {kern = new(SymmetricExponential)},
#                gaussian = {kern = new(Gaussian)},
#                powerlaw = {kern = new(PowerLaw)},
#                pareto3 = {kern = new(Pareto3)},
#                pareto2 = {kern = new(Pareto2)},
#                pareto1 = {kern = new(Pareto1)})
#     } else if ( # or that it refers to a valid hawkes kernel
#         !any(sapply(
#             paste0("Rcpp_", c("Exponential", "SymmetricExponential", "Gaussian", "PowerLaw", "Pareto3", "Pareto2", "Pareto1")),
#             function(class_) {is(kern, class_)}
#         ))
#     ) stop("'kern' must be a valid kernel.")
#
#     if (!is.null(binsize)) kern$binsize = binsize
#
#     # Periodogram
#     n <- length(counts)
#     dft <- fft(counts - mean(counts))
#     I <- Mod(dft)^2 / n
#     I <- I[-1]      # remove value corresponding to density in 0 (not well defined for centered processes)
#
#     # Whittle pseudo likelihood function (for optim)
#     nlopt_fn <- function(param) {
#         kern$param <- param
#         return( kern$whittle(I, trunc) )
#     }
#
#     if (is.null(opts))
#         opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "maxeval" = 1000)
#     else {
#         if (is.null(opts[["algorithm"]]))
#             opts <- c(opts, "algorithm" = "NLOPT_LN_NELDERMEAD")
#         if (is.null(opts[["maxeval"]]))
#             opts <- c(opts, "maxeval" = 1000)
#         if (is.null(opts[["xtol_rel"]]))
#             opts <- c(opts, "xtol_rel" = 1e-04)
#     }
#
#     # Sensible initialisation
#     if (is.null(init)) {
#         ymean = mean(counts)
#         # For PowerLaw
#         if (is(kern, "Rcpp_PowerLaw")) {
#             wmin = Inf
#             mu = .25
#             shape = 1
#             scale = 1
#             for (mu_ in c(.25, .5, .75)) {
#                 for (shape_ in 1:10) {
#                     for (scale_ in 1:10) {
#                         kern$param[1] = ymean * (1 - mu_)
#                         kern$param[2] = mu_
#                         kern$param[3] = shape_
#                         kern$param[4] = scale_
#                         whit = kern$whittle(I, trunc = trunc)
#                         if (whit < wmin) {
#                             mu = mu_
#                             shape = shape_
#                             scale = scale_
#                             wmin = whit
#                         }
#                     }
#                 }
#             }
#             x0 = c(ymean * (1 - mu), mu, shape, scale)
#         } else if (is(kern, "Rcpp_Exponential")) {
#             wmin = Inf
#             mu = .25
#             rate = 1
#             for (mu_ in c(.25, .5, .75)) {
#                 for (rate_ in 1:10) {
#                     kern$param[1] = ymean * (1 - mu_)
#                     kern$param[2] = mu_
#                     kern$param[3] = rate_
#                     whit = kern$whittle(I, trunc = trunc)
#                     if (whit < wmin) {
#                         mu = mu_
#                         rate = rate_
#                         wmin = whit
#                     }
#                 }
#             }
#             x0 = c(ymean * (1 - mu), mu, rate)
#         } else {
#             x0 = c(ymean * 2,
#                    .5,
#                    runif(length(kern$param)-2, 0, 10))
#         }
#     } else x0 = init
#
#     optargs = list(lb = rep(.0001, length(kern$param)),
#                    ub = c(Inf, .9999, rep(Inf, length(kern$param)-2)),
#                    opts = opts)
#     optargs = modifyList(optargs, list(...))
#     opt <- do.call(nloptr::nloptr, c(list(x0 = x0, eval_f = nlopt_fn), optargs))
#
#     # Create output object
#     kern$param = opt$solution
#
#     output = list(
#         par = opt$solution,
#         kernel = kern,
#         counts = counts,
#         binsize = kern$binsize,
#         opt = opt
#     )
#
#     class(output) = "HawkesModel"
#
#     return( output )
# }
