#' Fitting Hawkes processes from continuous data
#'
#' This function fits a Hawkes process to continuous data by minimizing the likelihood.
#'
#' @param model An object of class Model
#' @param data An object of class ContinuousData
#' @param opts (Optional) To be passed to `nloptr`
#' @param ... Additional arguments passed to `nloptr`
#'
#' @return
#' Returns the `optim` output
#' @export
#'
#' @examples
#' x = hawkes(100, fun = 1, repr = .5, family = "exp", rate = 1)
#' model = new(Exponential)
#' data = new(ContinuousData, x$p, x$T)
#' mle(model, data)
mle <- function(model, data, opts = NULL, ...) {

    model$attach(data)

    # Likelihood function (for nloptr)
    nlopt_fn <- function(param) {
        model$param <- param
        return( lapply(X = model$loglikngrad(), FUN = "-") )
    }

    if (is.null(opts))
        opts <- list("algorithm" = "NLOPT_LD_LBFGS")
    else if (is.null(opts[["algorithm"]]))
        opts <- c(opts, "algorithm" = "NLOPT_LD_LBFGS")

    optargs = list(lb = rep(.0001, length(model$param)),
                   ub = c(Inf, .9999, rep(Inf, length(model$param)-2)),
                   opts = opts)
    optargs = modifyList(optargs, list(...))
    opt <- do.call(nloptr::nloptr, c(list(x0 = model$param, eval_f = nlopt_fn), optargs))

    return( opt )
}
