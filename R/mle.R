#' Fitting Hawkes processes from continuous data
#'
#' This function fits a Hawkes process to continuous data by minimizing the likelihood
#' on the interval \eqn{[0,\mathrm{end}]}.
#'
#' @param events The locations of events (sorted in ascending order)
#' @param kern Either a string (partially) matching one of the kernels implemented (see Details), or an object of class Model
#' @param end The time until which the process is observed.
#' @param init (Optional) Initial values of the optimisation algorithm
#' @param opts (Optional) To be passed to `nloptr`
#' @param ... Additional arguments passed to `nloptr`
#'
#' @details
#' The maximum likelihood estimation procedure has only been implemented for the
#' exponential and the power law kernels.
#' For the exponential kernel, the likelihood is computed in \eqn{O(n)} complexity
#' (as described in details in T. Ozaki and Y. Ogata, “Maximum likelihood
#' estimation of Hawkes’ self-exciting point processes,” Ann. Inst. Stat. Math.,
#' vol. 31, no. 1, pp. 145–155, Dec. 1979).
#' For the power law kernel, the complexity is \eqn{O(n^2)}.
#'
#' @return
#' Returns a list containing the solution of the optimisation procedure, the object Model
#' with its parameters updated to the solution, and the output produced by `nloptr`.
#' @export
#'
#' @seealso [hawkes()] for the simulation of Hawkes processes,
#' [Model] for the abstract class, and [Exponential] for the specific
#' reproduction kernels.
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1,
#' # reproduction mean 0.5 and exponential fertility function with rate 2.
#' x = hawkes(100, fun = 1, repr = .5, family = "exp", rate = 1)
#' # Estimate the parameters from the arrival times of `x` using MLE
#' opt = mle(x$p, "Exponential", x$end)
#' opt$par                          # Estimated parameters
#' opt$model$ddloglik(x$p, x$end)     # Hessian matrix of the log-likelihood
mle <- function(events, kern, end, init = NULL, opts = NULL, ...) {

    # Check that the argument 'kern' is either a string that matches of the kernels implemented
    if (is.character(kern)) {
        kern = match.arg(tolower(kern),c("exponential", "powerlaw"))
        switch(kern,
               exponential = {model = new(Exponential)},
               powerlaw = {model = new(PowerLaw)})
    } else if ( # or that it refers to a valid hawkes model
        !any(sapply(
            paste0("Rcpp_", c("Exponential", "PowerLaw")),
            function(class_) {is(model, class_)}
        ))
    ) stop("'kern' must be a valid kernel.")

    # Likelihood function (for nloptr)
    nlopt_fn <- function(param) {
        model$param <- param
        return( lapply(X = model$loglikngrad(events, end), FUN = "-") )
    }

    if (is.null(opts))
        opts <- list("algorithm" = "NLOPT_LD_LBFGS")
    else {
        if (is.null(opts[["algorithm"]]))
            opts <- c(opts, "algorithm" = "NLOPT_LD_LBFGS")
        if (is.null(opts[["xtol_rel"]]))
            opts <- c(opts, "xtol_rel" = 1e-04)
    }

    if (is.null(init))
        x0 = c(runif(1, 0, 2),
               runif(1, 0, 1),
               runif(length(model$param)-2, 0, 2))
    else x0 = init

    optargs = list(lb = rep(.0001, length(model$param)),
                   ub = c(Inf, .9999, rep(Inf, length(model$param)-2)),
                   opts = opts)
    optargs = modifyList(optargs, list(...))
    opt <- do.call(nloptr::nloptr, c(list(x0 = x0, eval_f = nlopt_fn), optargs))

    # Create output object
    model$param = opt$solution

    output = list(
        par = opt$solution,
        model = model,
        events = events,
        end = end,
        opt = opt
    )

    class(output) = "HawkesModel"

    return( output )
}
